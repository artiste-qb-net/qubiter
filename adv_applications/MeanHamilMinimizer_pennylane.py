from adv_applications.MeanHamilMinimizer import *
from device_specific.Qubiter_to_PennyLane import *
import qml
import os
import importlib

class MeanHamilMinimizer_pennylane(MeanHamilMinimizer):

    def __init__(self, qasm_name, qnode_name, device, *args, **kwargs):

        MeanHamilMinimizer.__init__(self, *args, **kwargs)

        Qubiter_to_PennyLane(self.file_prefix, self.num_bits, self.hamil,
                             qasm_name, qnode_name, device, qasm_ftype='py')

        penny_file_path = self.file_prefix + '_' + qasm_name
        mod = importlib.import_module(penny_file_path)
        self.qnode_cost_fun = mod.getattr(mod, qnode_name)

    def find_min(self):
        x_val = self.init_x_val
        for k in range(self.num_iter):
            x_val = self.minimizer_fun.step(self.cost_fun, x_val)
        return x_val

    def cost_fun(self, x_val):
        cost = self.qnode_cost_fun(x_val)
        
        self.cur_x_val = x_val
        self.cur_cost = cost
        self.broadcast_cost_fun_call()
        self.iter_count += 1
            
        return cost