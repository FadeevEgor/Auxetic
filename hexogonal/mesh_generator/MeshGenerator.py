
from math import *
# N1, N2 --- количество элементов по горизонтали и вертикали соответственно
# l1 --- расстояние между узлами на основании элемента
# l2 --- расстояние между узлами по скошенной стороне
# fi --- угол между основание и скошенной сторонной в элементе
# x0, y0 ---  координаты левого нижнего узла сетки
n_long_rows = 0
length_of_row = 0
is_first_short = False
is_last_short = False
l1 = 0.0
l2 = 0.0
fi = 0.0
x0 = 0.0
y0 = 0.0


class Point2D:
    def __init__(self, x, y):
        self.x = x
        self.y = y


class Element:
    def __init__(self, indexes):
        self.index = indexes


class Mesh:
    def __init__(self):
        self.points = []
        self.elements = []


# чтение данных из файла
def read_input_data(file_name):
    global n_long_rows, length_of_row, l1, l2, fi, x0, y0, is_first_short, is_last_short
    in_file = open(file_name, 'r')
    input_data = in_file.read().split()
    in_file.close()
    n_long_rows = int(input_data[0])
    length_of_row = int(input_data[1])
    is_first_short = bool(int(input_data[2]))
    is_last_short = bool(int(input_data[3]))
    l1 = float(input_data[4])
    l2 = float(input_data[5])
    fi = grad_to_rad(float(input_data[6]))
    x0 = float(input_data[7])
    y0 = float(input_data[8])


def write_output_data(file_name, mesh):
    out_file = open(file_name, 'w')
    out_file.write(str(len(mesh.elements)) + ' ' + str(len(mesh.points)) + '\n')
    for element in mesh.elements:
        string = ''
        for index in element.index:
            string += str(index) + ' '
        out_file.write(string + '\n')
    for point in mesh.points:
        string = str(point.x) + ' ' + str(point.y)
        out_file.write(string + '\n')
    out_file.close()


def grad_to_rad(angle_in_grad):
    return angle_in_grad*pi/180.0


def row_of_points(x0, y0, dx, dy, n_elements):
    x = x0
    row = [Point2D(x0, y0)]
    for j in range(n_elements):
        x += dx
        row.append(Point2D(x, y0 + dy))
        x += dx
        row.append(Point2D(x, y0))
    return row


def long_row_of_elements(start_point_index, n_elements):
    elements = []
    for j in range(n_elements):
        left_down = start_point_index + 2 * j
        left_up = start_point_index + 2 * n_elements + 2 * j + 1
        indexes = [left_down, left_down + 1, left_down + 2,
                   left_up + 2, left_up + 1, left_up]
        elements.append(Element(indexes))
    return elements


def short_row_of_elements(start_point_index, n_elements):
    elements = []
    for j in range(n_elements - 1):
        left_down = start_point_index + 2*j
        left_up = start_point_index + 2*n_elements + 2*j + 1
        indexes = [left_down, left_down + 1, left_down + 2,
                   left_up + 2, left_up + 1, left_up]
        elements.append(Element(indexes))
    return elements


def first_row_of_elements(start_point_index, n_elements):
    elements = []
    for j in range(n_elements - 1):
        left_down = start_point_index + 2*j
        left_up = start_point_index + 2*n_elements + 2*j
        indexes = [left_down, left_down + 1, left_down + 2,
                   left_up + 2, left_up + 1, left_up]
        elements.append(Element(indexes))
    return elements


def last_row_of_elements(start_point_index, n_elements):
    elements = []
    for j in range(n_elements - 1):
        left_down = start_point_index + 2*j + 1
        left_up = start_point_index + 2*n_elements + 2*j + 1
        indexes = [left_down, left_down + 1, left_down + 2,
                   left_up + 2, left_up + 1, left_up]
        elements.append(Element(indexes))
    return elements


def get_mesh(mesh):
    dx = l2 * sin(fi)
    dy = l2 * cos(fi)
    delta_y = 2.0 * l1 - 2.0 * dy
    y = y0
    start_index = 0
    if is_first_short:
        mesh.points += row_of_points(x0 + dx, y0 - l1 + dy, dx, dy, length_of_row - 1)
        mesh.elements += first_row_of_elements(start_index, length_of_row)
        start_index += 2*length_of_row - 1

    for line in range(n_long_rows - 1):
        mesh.points += row_of_points(x0, y, dx, dy, length_of_row) + row_of_points(x0, y + l1, dx, -dy, length_of_row)
        y += delta_y
        mesh.elements += long_row_of_elements(start_index, length_of_row)
        start_index += 2*length_of_row + 1
        mesh.elements += short_row_of_elements(start_index + 1, length_of_row)
        start_index += 2*length_of_row + 1
    mesh.points += row_of_points(x0, y, dx, dy, length_of_row) + row_of_points(x0, y + l1, dx, -dy, length_of_row)
    mesh.elements += long_row_of_elements(start_index, length_of_row)
    start_index += 2 * length_of_row + 1

    if is_last_short:
        y += l1
        mesh.points += row_of_points(x0 + dx, y + l1 - dy, dx, -dy, length_of_row - 1)
        mesh.elements += last_row_of_elements(start_index, length_of_row)
        start_index += 2 * length_of_row - 1
    return mesh


read_input_data('in_data.txt')
mesh = get_mesh(Mesh())
output_file_name = 'mesh' + str(n_long_rows) + 'x' + str(length_of_row) + '.txt'
write_output_data(output_file_name, mesh)
