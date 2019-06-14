var_num_to_hist[1].append(None
)  # #1, line 1
for j20 in range(2):  # line 2
    var_num_to_hist[1].append(None
    )  # -my_fun1#1#2, line 3
    var_num_to_hist[2].append(None
    )  # -my_fun1#1#2, line 3
    fun_name_to_hist["my_fun1"].append(None
    )  # -my_fun1#1#2, line 3
    for j10 in range(4):  # line 4
        var_num_to_hist[1].append(None
        )  # -#1*.5, line 5
        var_num_to_hist[2].append(None
        )  # #2, line 5
        var_num_to_hist[3].append(None
        )  # -my_fun3#3, line 5
        fun_name_to_hist["my_fun3"].append(None
        )  # -my_fun3#3, line 5
    var_num_to_hist[1].append(None
    )  # my_fun1#1#2, line 7
    var_num_to_hist[2].append(None
    )  # my_fun1#1#2, line 7
    fun_name_to_hist["my_fun1"].append(None
    )  # my_fun1#1#2, line 7
var_num_to_hist[1].append(None
)  # #1*.3, line 9
var_num_to_hist[1].append(None
)  # my_fun#1, line 10
fun_name_to_hist["my_fun"].append(None
)  # my_fun#1, line 10

all_var_nums += [1, 2, 3]
all_fun_names += ['my_fun1', 'my_fun3', 'my_fun']
