import pandas as pan
import sys
if 'autograd.numpy' not in sys.modules:
    import numpy as np
import matplotlib.pyplot as plt
from qubiter.StateVec import *


class Plotter:
    """
    This class has no constructor. All its methods are static methods. It
    contains various utility methods for (1) creating pandas dataframes from
    state vectors and density matrices created by StateVec methods, and (2)
    plotting such dataframes. Dataframes with probability (resp., amplitude)
    entries are plotted as bar plots (resp., quiver plots). Casting things
    into dataframes is convenient because they are automatically displayed
    as html tables in jupyter notebooks.


    """
    @staticmethod
    def get_states(num_qbits, use_bin_labels=True, ZL=True):
        """
        Returns list of strings ['0', '1', ..] with 2^num_qbits entries if
        use_bin_labels=False, or the binary string equivalents if
        use_bin_labels=True.

        Parameters
        ----------
        num_qbits : int
        use_bin_labels : bool
        ZL : bool
            If True, appends "(ZL)" to state labels. Else appends "(ZF)"

        Returns
        -------
        list[str]

        """
        str1 = ""
        if num_qbits > 1:
            if ZL:
                str1 = "ZL"
            else:
                str1 = "ZF"
        if not use_bin_labels:
            states = [str(x) for x in range(0, 1 << num_qbits)]
        else:
            states = [np.binary_repr(x, width=num_qbits) + str1
                      for x in range(0, 1 << num_qbits)]
        return states

    @staticmethod
    def get_st_vec_df(trad_st_vec, use_bin_labels=True):
        """
        Admits as input trad_st_vec which is a complex numpy array with
        shape (dim,) where dim= 2^num_qbits and ZL convention is assumed.
        Returns a dataframe with a single column with 2^num_bit rows. The
        rows are labelled by 0, 1, 2, ... or their binary equivalents,
        depending on the bool value of use_bin_labels.

        Parameters
        ----------
        trad_st_vec : np.ndarray
        use_bin_labels : bool

        Returns
        -------
        pan.DataFrame

        """
        num_qbits = trad_st_vec.shape[0].bit_length() - 1
        states = Plotter.get_states(num_qbits, use_bin_labels)
        return pan.DataFrame(trad_st_vec, index=states)

    @staticmethod
    def get_den_mat_df(num_qbits, den_mat, use_bin_labels=True):
        """
        Admits as input a numpy array den_mat (a density matrix) of shape (
        dim, dim), where dim= 2^num_qbits and ZL convention is assumed.
        Returns square dataframe with entries of den_mat. The rows and
        columns are labelled by 0, 1, 2, ... or their binary equivalents,
        depending on the bool value of use_bin_labels.

        Parameters
        ----------
        num_qbits : int
        den_mat : np.ndarray
        use_bin_labels : bool

        Returns
        -------
        pan.DataFrame

        """
        assert 1 << num_qbits == den_mat.shape[0]
        states = Plotter.get_states(num_qbits, use_bin_labels)
        return pan.DataFrame(den_mat, columns=states, index=states)

    @staticmethod
    def get_pd_df(num_qbits, pd, use_bin_labels=True):
        """
        Admits as input a 1-dim numpy array pd (probability distribution)
        with shape ( 2^num_qbits, ) and ZL convention is assumed. Returns a
        dataframe with its entries. The rows are labelled by 0, 1, 2,
        ... or their binary equivalents, depending on the bool value of
        use_bin_labels.

        Parameters
        ----------
        num_qbits : int
        pd : np.ndarray
        use_bin_labels : bool

        Returns
        -------
        pan.DataFrame

        """
        # print(pd)
        assert 1 << num_qbits == pd.shape[0]
        states = Plotter.get_states(num_qbits, use_bin_labels)
        return pan.DataFrame(pd, index=states)

    @staticmethod
    def get_bit_probs_df(bit_probs):
        """
        Admits as input a list bit_probs with num_qbits items consisting of
        pairs of two floats. Returns a dataframe with two columns and
        num_qbits rows. The rows are labelled by 0, 1, 2, ...
      
        Parameters
        ----------
        bit_probs : list[tuple(float, float)]

        Returns
        -------
        pan.DataFrame

        """
        p0_col, p1_col = zip(*bit_probs)
        # print('p0_col', p0_col)
        # print('p1_col', p1_col)
        return pan.DataFrame({'Prob(0)': p0_col,
                              'Prob(1)': p1_col})
    
    @staticmethod
    def plot_probs_col(titles, probs_col_df_list):
        """
        Admits as input a list of dataframes called 'probs_col_df_list'. The
        names of the dataframes are given by the list of strings 'titles'.
        The dataframes all have a single column, and possibly different
        numbers of rows. All entries of all dataframes are floats between 0
        and 1. This function plots each dataframe in the list as a barplot.

        Parameters
        ----------
        titles : list[str]
        probs_col_df_list : list[pan.DataFrame]

        Returns
        -------
        None

        """

        assert len(titles) == len(probs_col_df_list)
        num_titles = len(titles)

        def single_pd(ax, title, pd_df):
            y_pos = np.arange(len(pd_df.index)) + .5
            plt.sca(ax)
            plt.yticks(y_pos, pd_df.index)
            ax.invert_yaxis()

            ax.set_xticks([0, .25, .5, .75, 1])
            ax.set_xlim(0, 1)

            for row in range(len(y_pos)):
                val = pd_df.values[row]
                if isinstance(val, np.ndarray):
                    val = val[0]
                ax.text(val, y_pos[row], '{:.3f}'.format(val))

            ax.grid(True)
            ax.set_title(title)
            # new version of python/matplotlib has bug here.
            # The following used to work but no longer does.
            # ax.barh(y_pos, pd_df.values, align='center')
            # work around
            for b in range(len(y_pos)):
                ax.barh(y_pos[b], pd_df.values[b],
                        align='center', color='blue')
        plt.close('all')
        fig, ax_list = plt.subplots(nrows=num_titles, ncols=1)
        # print("***", ax_list[0])
        if num_titles == 1:
            ax_list = [ax_list]
        for k, tit in enumerate(titles):
            single_pd(ax_list[k], tit, probs_col_df_list[k])
        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_phasors(titles, st_vec_df_list=None, 
                     den_mat_df_list=None):
        """
        Admits as input a list of dataframes called 'st_vec_df_list' or one
        called 'den_mat_df_list' but not both. Exactly one of these lists
        must be None. The names of the dataframes are given by the list of
        strings 'titles'. Let dim = 2^num_qbits. The dataframes in
        st_vec_df_list have one column with dim rows, whereas those for
        den_mat_df_list are square with dim rows and columns. This function
        plots each dataframe in the list as a quiver plot (phasor->arrow)
        with dim rows and either one column or dim columns.

        Parameters
        ----------
        titles : list[str]
        st_vec_df_list : list[pan.DataFrame]|None
        den_mat_df_list : list[pan.DataFrame]|None

        Returns
        -------
        None

        """
        flag = 0
        if st_vec_df_list is not None:
            flag += 1
        if den_mat_df_list is not None:
            flag += 1
        assert flag == 1
        if den_mat_df_list:
            case = 'dm'
            df_list = den_mat_df_list
        else:
            case = 'vec'
            df_list = st_vec_df_list
        num_titles = len(titles)
        
        def single_ax(ax, title, df):
            states = df.index
            num_sts = len(states)
            y = np.linspace(0, num_sts-1, num_sts)
            if case == 'vec':
                assert len(df.columns) == 1
                x = [0]
                x_num_sts = 1
            else:
                assert len(df.columns) == num_sts
                x = y
                x_num_sts = num_sts
            xx, yy = np.meshgrid(x, y)
            # print(xx)
            # print(yy)

            ax.set_xlim(-1, x_num_sts)
            ax.set_xticks(np.arange(0, x_num_sts))
            ax.xaxis.tick_top()

            ax.set_ylim(-1, num_sts)
            ax.set_yticks(np.arange(0, num_sts))
            ax.invert_yaxis()

            ax.set_aspect('equal', adjustable='box')
            ax.set_title(title, y=1.2)

            for st, nom in enumerate(states):
                ax.annotate(nom, xy=(x_num_sts+.25, st),
                            annotation_clip=False)
            max_mag = np.max(np.absolute(df.values))
            q = ax.quiver(xx, yy, df.values.real,
                df.values.imag, scale=max_mag, units='x')

            plt.quiverkey(q, 0, -2, max_mag,
                '= {:.2e}'.format(max_mag), labelpos='E', coordinates='data')

        plt.close('all')
        fig, ax_list = plt.subplots(nrows=num_titles, ncols=1)
        if num_titles == 1:
            ax_list = [ax_list]
        for k, tit in enumerate(titles):
            single_ax(ax_list[k], tit, df_list[k])
        plt.tight_layout(pad=2)
        plt.show()

    @staticmethod
    def plot_counts(state_name_to_counts):
        """
        Plots an OrderedDict[str, int] called state_name_to_counts as a bar
        graph. Divides counts by their sum before plotting.

        Parameters
        ----------
        state_name_to_counts : OrderedDict[str, int]

        Returns
        -------
        None

        """
        dictio = state_name_to_counts
        data = np.array(list(dictio.values()), dtype=float)
        data /= np.sum(data)
        df = pan.DataFrame(data, index=dictio.keys())
        Plotter.plot_probs_col(['count prob'], [df])


if __name__ == "__main__":

    def main():
        num_qbits = 3
        st_vec0 = StateVec(num_qbits,
            arr=StateVec.get_random_st_vec(num_qbits).arr)
        st_vec1 = StateVec(num_qbits,
            arr=StateVec.get_random_st_vec(num_qbits).arr)
        st_vec_dict = {'br0': st_vec0,
                       'br1': st_vec1,
                       'br3': None}

        trad_st_vec = st_vec0.get_traditional_st_vec
        den_mat = StateVec.get_den_mat(num_qbits, st_vec_dict)
        # print("den_mat", den_mat)
        st_vec_pd = st_vec0.get_pd()
        den_mat_pd = StateVec.get_den_mat_pd(den_mat)
        bit_probs_vec = StateVec.get_bit_probs(num_qbits, st_vec_pd)
        bit_probs_dm = StateVec.get_bit_probs(num_qbits, den_mat_pd)

        st_vec_df = Plotter.get_st_vec_df(st_vec0.get_traditional_st_vec)
        den_mat_df = Plotter.get_den_mat_df(num_qbits, den_mat)
        # print("den_mat_df", den_mat_df)
        st_vec_pd_df = Plotter.get_pd_df(num_qbits, st_vec_pd)
        den_mat_pd_df = Plotter.get_pd_df(num_qbits, den_mat_pd)
        bit_probs_df1 = Plotter.get_bit_probs_df(bit_probs_vec)
        bit_probs_df2 = Plotter.get_bit_probs_df(bit_probs_dm)

        Plotter.plot_probs_col(['st_vec_pd'], [st_vec_pd_df])
        # print(bit_probs_df1)
        Plotter.plot_probs_col(
            ['bit_probs, Prob(0)'], [bit_probs_df1["Prob(0)"]])
        Plotter.plot_phasors(['st_vec'], st_vec_df_list=[st_vec_df])
        Plotter.plot_phasors(['den_mat'], den_mat_df_list=[den_mat_df])
    main()
