"""
NOTA    |   |   |   |   |   |
NOTA    |  |'''''''''|  |   |
NOTA    |  |Message  |  |   |
NOTA    |  |.........|  |   |
NOTA    |   |   |   |   |   |

NOTA    |   @---+---@---@   |
"""

def NOTA_box_str(num_bits, start_bit, end_bit, message):
    assert start_bit > end_bit
    assert num_bits > start_bit
    assert end_bit > -1
    str = "NOTA\t|"
    cur_bit = 0
    while cur_bit > start_bit:
        str += "   |"
    str += "  |'"
