from qubiter.SEO_writer import *
import itertools as it
import pprint as pp
import collections as col


class StairsCkt_writer(SEO_writer):
    """
    This class is a subclass of class SEO_writer and it writes a "Stairs
    Circuit". For example, this is what the Picture file of a Stairs Circuit 
    looks like for num_qbits = 3
    
    U   |   |
    O---U   |  
    @---U   |   
    O---O---U   
    O---@---U   
    @---O---U  
    @---@---U 
       
    Here, U is a general U(2) matrix with 4 parameters, all of which can be 
    made into placeholder variables. If each U is represented by a node and 
    the controls of each U represent its parents, then this quantum circuit 
    can be represented by a fully connected Quantum Bayesian Network (QB 
    net). (See my >10 year old blog called "Quantum Bayesian Networks" for 
    more info than you would ever want to know about QB nets). 
    
    This class can also be asked to construct a QB net that is **not** fully 
    connected, by limiting the number of controls for a given U to fewer 
    than all the ones to its left. For example, suppose that in the 
    num_qbits=3 case, we restrict the parents of the U in the last step to
    just one, instead of the 2 parents that it has in the fully connected 
    case. Then we get 

    
    U   |   |
    O---U   |  
    @---U   |   
    O---+---U   
    @---+---U  
       
    or 
    
    U   |   |
    O---U   |  
    @---U   |   
    |   O---U   
    |   @---U     
    
    The constructor of this class has as input an ordered dictionary called 
    gate_str_to_rads_list. This dictionary gives for each gate in the 
    quantum circuit, a gate string gate_str that specifies the gate. 
    gate_str_to_rads_list maps gate_str to a list of 4 floats (or 
    placeholder variables for those floats) for the 4 parameters of the U 
    matrix. For example, here are possible values for gate_str_to_rads_list 
    for the num_qbits=3 fully connected qb net
    
    with every rads_list item filled with the same constant .3
    {'prior': [0.3, 0.3, 0.3, 0.3],
     '2F': [0.3, 0.3, 0.3, 0.3],
     '2T': [0.3, 0.3, 0.3, 0.3],
     '2F1F': [0.3, 0.3, 0.3, 0.3],
     '2F1T': [0.3, 0.3, 0.3, 0.3],
     '2T1F': [0.3, 0.3, 0.3, 0.3],
     '2T1T': [0.3, 0.3, 0.3, 0.3]}
     
    with every rads_list item filled by a random number between 0 and 2pi
    {'prior': [0.46731839721496604,
               0.012285135138256131,
               0.20001353832948487,
               0.36694428209569985],
     '2F': [4.1968011007222898,
            5.1978252498063808,
            4.8063090848060321,
            4.2509081392354409],
     '2T': [4.3359074640905213,
            2.0749617893052315,
            4.555666727197961,
            5.3092010293653802],
     '2F1F': [0.99177045463186475,
              3.3344615340103325,
              2.1441702948866386,
              2.4603764283165521],
     '2F1T': [4.0909522483111145,
              2.0714182784661888,
              5.4034187072431923,
              6.0856723571386766],
     '2T1F': [4.0000452017061194,
              3.7193341571216658,
              3.381322125034953,
              5.4492142181489802],
     '2T1T': [6.2597553541046853,
              0.077807529496169509,
              3.7389318319862217,
              6.2233264819972307]}
              
    with every rads_list item filled by a unique placeholder variable string
    {'prior': ['#50', '#51', '#52', '#53'],
     '2F': ['#500', '#501', '#502', '#503'],
     '2T': ['#510', '#511', '#512', '#513'],
     '2F1F': ['#5000', '#5001', '#5002', '#5003'],
     '2F1T': ['#5010', '#5011', '#5012', '#5013'],
     '2T1F': ['#5100', '#5101', '#5102', '#5103'],
     '2T1T': ['#5110', '#5111', '#5112', '#5113']}

    This is what gate_str_to_rads_list looks like in the num_qbits=3 case,
    when the last U has only one parent (qbit 2) instead of two parents (
    qbits 1 and 2): 
    
    {'prior': ['#50', '#51', '#52', '#53'],
     '2F': ['#500', '#501', '#502', '#503'],
     '2T': ['#510', '#511', '#512', '#513'],
     '2F1_': ['#5050', '#5051', '#5052', '#5053'],
     '2T1_': ['#5150', '#5151', '#5152', '#5153']}
     
    Note that all placeholder strings begin with '#5' to insure that once
    the hash character is removed, the remaining number doesn't start with
    '0'. Note that characters '_' and '5' represent bits whose values are
    unspecified.

    Attributes
    ----------
    gate_str_to_rads_list : OrderedDict[str, list[float|str]]

    """
    
    def __init__(self,  gate_str_to_rads_list, file_prefix, emb, **kwargs):
        """
        Constructor

        This constructor writes English and Picture files but it doesn't
        close those files after writing them. You must do that yourself
        using close_files().


        Parameters
        ----------
        gate_str_to_rads_list : dict[str, list[float|str]]
        file_prefix : str
            file prefix for English and Picture files written by this class
        emb : CktEmbedder
        kwargs : dict
            key-word arguments of SEO_writer

        Returns
        -------

        """
        SEO_writer.__init__(self, file_prefix, emb, **kwargs)
        self.gate_str_to_rads_list = gate_str_to_rads_list
        self.write()

    @staticmethod
    def get_gate_str_to_rads_list(num_qbits, fill_type, rads_const=None,
                                  u2_bit_to_higher_bits=None):
        """
        This method returns a gate_str_to_rads_list constructed according to
        the specs given by its arguments.

        fill_type is a string in ['const', 'rand', '#int'] The 3 types of
        fill_type have already been illustrated in the class docstring. If
        the fill_type is 'const', then the method expects a float for
        rads_const.

        u2_bit_to_higher_bits is used to restrict the controls of each U.
        For example, for num_qbits=3,

        u2_bit_to_higher_bits = {0: [1, 2], 1: [2], 2: []}

        specifies a fully connected qb net, whereas

        u2_bit_to_higher_bits = {0: [2], 1: [2], 2: []}

        means qubit 0 has qubit 2 but not 1 as parent.

        Parameters
        ----------
        num_qbits : int
        fill_type : str
            either 'const', 'rand' or '#int'
        rads_const : float | None
        u2_bit_to_higher_bits : dict[int, list[int]]

        Returns
        -------
        OrderedDict[str, list[float]]

        """
        # each "step" may have several gate_str of same length
        const_list = [rads_const for k in range(4)]
        gate_str_to_rads_list = col.OrderedDict()
        if fill_type == 'const':
            assert rads_const is not None
            gate_str_to_rads_list['prior'] = const_list
        elif fill_type == 'rand':
            rand_list = list(np.random.random_sample((4,)))
            gate_str_to_rads_list['prior'] = rand_list
        elif fill_type == '#int':
            # all placeholder variables will start
            # with 5 to avoid starting with 0
            gate_str_to_rads_list['prior'] = ['#50', '#51', '#52', '#53']
        else:
            assert False, 'unsupported fill type'

        pair = ['F', 'T']
        singlet = ['_']
        for tup_len in range(1, num_qbits):
            u2_pos = num_qbits - tup_len - 1
            pa_range = range(u2_pos+1, num_qbits)
            parent_to_list = {k: pair for k in pa_range}
            if u2_bit_to_higher_bits:
                parent_to_list = {k: singlet for k in pa_range}
                for pa_bit in u2_bit_to_higher_bits[u2_pos]:
                    assert u2_pos < pa_bit < num_qbits
                    parent_to_list[pa_bit] = pair
            list_of_lists = [parent_to_list[k] for k in reversed(pa_range)]
            # print("mmmnnnnnn", list_of_lists)
            for tuple_of_FTs in it.product(*list_of_lists):
                s = ''
                for k in range(tup_len):
                    s += str(num_qbits - 1 - k) + tuple_of_FTs[k]
                if fill_type == 'const':
                    gate_str_to_rads_list[s] = const_list
                elif fill_type == 'rand':
                    rand_list = list(2*np.pi*np.random.random_sample((4,)))
                    gate_str_to_rads_list[s] = rand_list
                elif fill_type == '#int':
                    hash_str = '#5'
                    for k in range(tup_len):
                        x = tuple_of_FTs[k]
                        if x == '_':
                            hash_str += '5'
                        elif x == 'F':
                            hash_str += '0'
                        elif x == 'T':
                            hash_str += '1'
                    gate_str_to_rads_list[s] =\
                        [hash_str + str(k) for k in range(4)]
                else:
                    assert False
        return gate_str_to_rads_list

    @staticmethod
    def get_all_var_nums(gate_str_to_rads_list):
        """
        This method scans each rads_list of gate_str_to_rads_list for items
        of the type '#x' where x is an int. Every int x is added to a list
        all_var_nums which is returned.

        Parameters
        ----------
        gate_str_to_rads_list : OrderedDict[str, list[float|str]]

        Returns
        -------
        list[int]

        """
        all_var_nums = []
        for rads_list in gate_str_to_rads_list.values():
            for rads in rads_list:
                if isinstance(rads, str) and rads[0] == '#':
                    all_var_nums.append(int(rads[1:]))
        return all_var_nums

    @staticmethod
    def get_var_num_to_rads(gate_str_to_rads_list, fill_type,
                            rads_const=None):
        """
        This method returns a dict var_num_to_rads obtained as follows. The
        rads lists in gate_str_to_rads_list are scanned for items of the
        type '#x' where x is an int. Then x is mapped to a float, either the
        constant rads_const if fill_type is 'const', or a random number if
        fill type is 'rand'.

        Parameters
        ----------
        gate_str_to_rads_list : OrderedDict[str, list[float|str]]
        fill_type : str
        rads_const : float | None

        Returns
        -------
        dict[int, float]

        """
        var_num_to_rads = {}
        for rads_list in gate_str_to_rads_list.values():
            for rads in rads_list:
                if isinstance(rads, str) and rads[0] == '#':
                    if fill_type == 'const':
                        assert rads_const is not None
                        var_num_to_rads[int(rads[1:])] = rads_const
                    elif fill_type == 'rand':
                        var_num_to_rads[int(rads[1:])] =\
                            2*np.pi*np.random.random()
                    else:
                        assert False, 'unsupported fill type'
        return var_num_to_rads

    @staticmethod
    def make_array_from_gate_str_to_rads_list(gate_str_to_rads_list):
        """
        This method returns a numpy array which is constructed from
        gate_str_to_rads_list by vertically stacking (with np.vstack()) its
        rads_lists

        Parameters
        ----------
        gate_str_to_rads_list : dict[str, list[float|str]]

        Returns
        -------
        np.ndarray

        """
        return np.vstack([np.array(rads_list) for rads_list in
                          gate_str_to_rads_list.values()])

    def get_u2_pos(self, gate_str):
        """
        Given a well formed gate_str (one of the keys of
        gate_str_to_rads_list), this method returns the bit position of the
        U(2) matrix.

        Parameters
        ----------
        gate_str : str

        Returns
        -------
        int

        """
        num_qbits = self.emb.num_qbits_bef
        if gate_str != 'prior':
            u2_pos = num_qbits - len(gate_str) // 2 - 1
        else:
            u2_pos = num_qbits-1
        return u2_pos

    @staticmethod
    def get_controls_from_gate_str(num_qbits, gate_str):
        """
        This method returns an object of class Controls, constructed from
        info in the input `gate_str` (a well formed key of
        gate_str_to_rads_lis)


        Parameters
        ----------
        num_qbits : int
        gate_str : str

        Returns
        -------
        Controls

        """
        trols = Controls(num_qbits)
        if gate_str != 'prior':
            for k in range(len(gate_str)//2):
                trol_pos = int(gate_str[2 * k])
                trol_kind = gate_str[2 * k + 1]
                # allow for possibility that trol_kind = '_' (no trol)
                if trol_kind == 'T':
                    trols.bit_pos_to_kind[trol_pos] = True
                elif trol_kind == 'F':
                    trols.bit_pos_to_kind[trol_pos] = False
        trols.refresh_lists()
        return trols
        
    def write(self):
        """
        This method writes English and Picture files for a Stairs Circuit.

        Returns
        -------

        """
        num_qbits = self.emb.num_qbits_bef
        for gate_str, rads_list in self.gate_str_to_rads_list.items():
            num_qbits = self.emb.num_qbits_bef
            trols = StairsCkt_writer.get_controls_from_gate_str(
                num_qbits, gate_str)
            u2_pos = self.get_u2_pos(gate_str)
            self.write_controlled_one_bit_gate(u2_pos, trols,
                    OneBitGates.u2, rads_list)


if __name__ == "__main__":
    def main():
        num_qbits = 3
        for fill_type in ['const', 'rand', '#int']:
            di = StairsCkt_writer.get_gate_str_to_rads_list(
                num_qbits, fill_type, rads_const=.3)
            pp.pprint(di)
         
        u2_bit_to_higher_bits = {0: [2], 1: [2], 2: []}
        di = StairsCkt_writer.get_gate_str_to_rads_list(
                num_qbits, "#int", u2_bit_to_higher_bits=u2_bit_to_higher_bits)
        pp.pprint(di)

        vn_to_r = StairsCkt_writer.get_var_num_to_rads(di,
                                                       fill_type='const',
                                                       rads_const=.3)
        pp.pprint(vn_to_r)

        arr = StairsCkt_writer.make_array_from_gate_str_to_rads_list(di)
        print("arr=\n", arr)

        num_qbits = 4
        gate_str_to_rads_list = StairsCkt_writer.get_gate_str_to_rads_list(
            num_qbits, '#int')
        file_prefix = 'stairs_writer_test'
        emb = CktEmbedder(num_qbits, num_qbits)

        wr = StairsCkt_writer(gate_str_to_rads_list, file_prefix, emb)
        wr.close_files()
        wr.print_eng_file()
        wr.print_pic_file()

    main()
