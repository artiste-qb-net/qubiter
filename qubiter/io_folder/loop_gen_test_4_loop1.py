var_num_to_hist[1].append(3
)  # #1, line 1
for j20 in range(2):  # line 2
    fun_name_to_hist["my_fun1"].append(
        lambda x1, x2, j20=j20: x1 + x2 + j20
    )  # -my_fun1#1#2, line 3
    var_num_to_hist[1].append(j20*.3
    )  # -my_fun1#1#2, line 3
    var_num_to_hist[2].append(.45
    )  # -my_fun1#1#2, line 3
    for j10 in range(4):  # line 4
        var_num_to_hist[1].append(j10 + j20*.2
        )  # -#1*.5, line 5
        var_num_to_hist[2].append(None
        )  # #2, line 5
        fun_name_to_hist["my_fun3"].append(
            lambda x3, j10=j10, j20=j20: -x3*.5 + j10 + j20
        )  # -my_fun3#3, line 5
        var_num_to_hist[3].append(2*j10
        )  # -my_fun3#3, line 5
    fun_name_to_hist["my_fun1"].append(None
    )  # my_fun1#1#2, line 7
    var_num_to_hist[1].append(None
    )  # my_fun1#1#2, line 7
    var_num_to_hist[2].append(None
    )  # my_fun1#1#2, line 7
var_num_to_hist[1].append(None
)  # #1*.3, line 9
fun_name_to_hist["my_fun"].append(
    lambda x1: x1*.7
)  # my_fun#1, line 10
var_num_to_hist[1].append(None
)  # my_fun#1, line 10

all_var_nums += [1, 2, 3]
all_fun_names += ['my_fun1', 'my_fun3', 'my_fun']

# import pprint as pp
# print("---------------")
# print("var_num_to_hist")
# pp.pprint(var_num_to_hist)
# print()
# print("fun_name_to_hist")
# pp.pprint(fun_name_to_hist)