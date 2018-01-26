
ibmqx2_edges = (
    (0, 2),
    (0, 1),
    (1, 2),
    (3, 2),
    (3, 4),
    (4, 2))  # 6 edges

ibmqx4_edges = (
    (1, 0),
    (2, 1),
    (2, 0),
    (2, 4),
    (3, 2),
    (3, 4))  # 6 edges

ibmqx5_edges = (
    (1, 2),
    (1, 0),
    (2, 3),
    (3, 4),
    (3, 14),
    (5, 4),
    (6, 5),
    (6, 7),
    (6, 11),
    (7, 10),
    (8, 7),
    (9, 8),
    (9, 10),
    (11, 10),
    (12, 5),
    (12, 13),
    (12, 11),
    (13, 14),
    (13, 4),
    (15, 0),
    (15, 2),
    (15, 14))  # 22 edges

# Simulation only --- all CNOT connections are allowed
def simulator_edges(num_bits):
    edges = ()
    for i in range(num_bits):
        for j in range(num_bits):
            if i != j:
                edges += ((i, j),)
    return edges


if __name__ == '__main__':
    print(simulator_edges(4))


