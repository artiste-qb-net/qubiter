import pyclbr
import os
import pprint as pp
import networkx as nx
import graphviz as gv
from shutil import copyfile

"""
The purpose of this script is to generate a file called classgraph.pdf 
with the class graph for qubiter. The script also generates a file called 
classgraph_orphans.html that lists the orphan classes. Orphan classes (i.e., 
classes with no parents or children) don't show up in classgraph.pdf. 
"""

# cl_to_data = pyclbr.readmodule("SEO_reader")
# print(name_to_data)

# create cl_to_parents: dict[str, list[str]]
cl_to_parents = {}

dir_whitelist = [
    "./",
    'adv_applications',
    'device_specific',
    'latex_tools',
    'quantum_CSD_compiler'
    ]
file_blacklist = [
    '__init__.py',
    'run_all_nb.py',
    'run_all_py.py',
    'classgraph.py'
]
for dir_name in dir_whitelist:
    for fname in os.listdir(dir_name):
        if fname[-3:] == '.py' and fname not in file_blacklist:
            # print("mmmnnnn", dir_name, fname)
            cl_to_data = pyclbr.readmodule(fname[:-3],
                        [os.path.abspath(dir_name)])
            for cl_name, cl_desc in cl_to_data.items():
                if os.path.basename(cl_desc.file) != fname:
                    continue
                pa_list = [x.name for x in cl_desc.super]
                cl_to_parents[cl_name] = pa_list
pp.pprint(cl_to_parents)

# create nx graph and orphans list from cl_to_parents
nx_graph = nx.DiGraph()
roots = []
for nd, pa_list in cl_to_parents.items():
    for pa in pa_list:
        nx_graph.add_edge(pa, nd)
    if not pa_list:
        roots.append(nd)
orphans = []
for nd in roots:
    if nd not in nx_graph:
        orphans.append(nd)
orphans.sort()

# write dot file from nx_graph
nom = 'classgraph'
dot_name = nom + '.dot'
nx.nx_pydot.write_dot(nx_graph, dot_name)

# write pdf file from dot file
with open(dot_name, 'r') as fi:
    src = gv.Source(fi.read())
src.render(dot_name, view=False, format='pdf')
pdf_name = nom + '.pdf'
copyfile(dot_name + '.pdf', pdf_name)
os.remove(dot_name + '.pdf')

# write txt file listing orphans (not shown in classgraph)
orp_name = 'classgraph_orphans.html'
with open(orp_name, 'w') as fi:
    fi.write('<html><body>\n')
    fi.write('<ol>\n')
    for nd in orphans:
        fi.write('<li>' + nd + '\n')
    fi.write('</ol>\n')
    fi.write('</body></html>')





