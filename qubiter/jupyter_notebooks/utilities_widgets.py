import ipywidgets as wid
from IPython.display import display, clear_output
import numpy as np
from qubiter.PlaceholderManager import *
from qubiter.SEO_simulator import *
from qubiter.StateVec import *


def run_sim_gui(file_prefix, num_bits, all_var_nums, fun_name_to_fun=None,
                slider_max_degs=360*3, append_new=False,
                sty_fin_desc='ALL'):
    """
    Generates and runs a widgets gui (graphical user interface). The gui has
    a button labelled `Run` that creates an object of SEO_simulator. The gui
    contains a slider for each placeholder variable (parameter) of a circuit
    that has been created a priori by a SEO_writer using the name
    `file_prefix` and number of qubits `num_bits`.

    If a parameter is labelled `#1`, then the slider value degs_1 times
    pi/180 is substituted for `#1`

    Parameters
    ----------
    file_prefix : str
    num_bits : int
    all_var_nums : list[int]
        all the placeholder variable numbers. If the circuit has exactly
        two placeholder variables #5 and #7, then all_var_nums=[5, 7]
    fun_name_to_fun : dict[str, function]
        dict mapping every functional placeholder name to its function
    slider_max_degs : float
        maximum (in degrees) of sliders (same for all of them)
    append_new : bool
        If True, printout for current run will be appended to end. If False,
        the previous printout will be erased before printing for current run.
    sty_fin_desc : str
        The style used in the description of the final state vector.
        Argument of StateVec.get_style_dict().

    Returns
    -------
    None

    """
    var_num_to_rads = {}
    hbox_comps = []
    slider_list = []
    for num in all_var_nums:
        slider = wid.FloatSlider(
            min=0.0,
            max=slider_max_degs,
            step=10,
            orientation='vertical')
        slider_list.append(slider)
        text = wid.FloatText(layout=wid.Layout(width='70px'))
        label = wid.Label('degs_' + str(num))
        wid.jslink((slider, 'value'), (text, 'value'))
        hbox_comps.append(wid.VBox([slider, text, label]))
        var_num_to_rads[num] = text.value * np.pi / 180

    rbut = wid.Button(description='Run')
    rbut.layout.width = '110px'
    rbut.button_style = 'danger'

    ebut = wid.Button(description='Erase Output')
    ebut.layout.width = '110px'
    ebut.button_style = 'danger'

    hbox_comps.append(wid.VBox([rbut, ebut]))
    hbox = wid.HBox(hbox_comps)
    display(hbox)

    def slider_do(num1, change):
        if change['new'] != change['old']:
            var_num_to_rads[num1] = change['new'] * np.pi / 180
        else:
            var_num_to_rads[num1] = change['old'] * np.pi / 180

    for num1, slider in zip(all_var_nums, slider_list):
        fun = (lambda x, num2=num1: slider_do(num2, x))
        slider.observe(fun, names='value')

    def rbut_event_handler(btn):
        vman = PlaceholderManager(var_num_to_rads=var_num_to_rads,
                                  fun_name_to_fun=fun_name_to_fun)
        if not append_new:
            clear_output()
            display(hbox)
        for num1 in all_var_nums:
            rads = var_num_to_rads[num1]
            print('degs_' + str(num1) + '(=rads)',
                  rads*180/np.pi, '(',  rads, ')')
        sim = SEO_simulator(file_prefix, num_bits, vars_manager=vman)
        print('\n-----------------------------beginning final results')
        StateVec.describe_st_vec_dict(sim.cur_st_vec_dict,
                                      **StateVec.get_style_dict(sty_fin_desc))
        print('-----------------------------ending final results\n')
    rbut.on_click(rbut_event_handler)

    def ebut_event_handler(btn):
        clear_output()
        display(hbox)
    ebut.on_click(ebut_event_handler)
