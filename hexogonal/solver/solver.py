# TODO:
# refactoring
import numpy as np
from scipy import optimize as opt


# parameters of functional
w_a = 0.3
w_l = 1.0
# points to calculate Poison' ratio
bottom = (28, 41)
top = (36, 49)
right = 6
left = 71


class Point2D:
    def __init__(self, x=0.0, y=0.0):
        self.x, self.y = x, y

    def set_coordinates(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return str(self.x) + " " + str(self.y)

    def __add__(self, other):
        return Point2D(self.x + other.x, self.y + other.y)


class Rib:
    def __init__(self, indices, length):
        self.indices = [index for index in indices]
        self.l0 = length


class Element:
    def __init__(self, indices):
        self.indices = [index for index in indices]

    def set_indices(self, indices):
        self.indices = [index for index in indices]

    def angles(self):
        result = np.zeros(dtype=np.int32, shape=[4, 3])
        result[0, :] = [self.indices[0], self.indices[-1], self.indices[1]]
        result[1, :] = [self.indices[2], self.indices[1], self.indices[3]]
        result[2, :] = [self.indices[3], self.indices[2], self.indices[4]]
        result[3, :] = [self.indices[5], self.indices[4], self.indices[0]]
        return result

    def ribs(self):
        result = np.zeros(dtype=np.int32, shape=[6, 2])
        result[0, :] = [self.indices[0], self.indices[1]]
        result[1, :] = [self.indices[1], self.indices[2]]
        result[2, :] = [self.indices[2], self.indices[3]]
        result[3, :] = [self.indices[3], self.indices[4]]
        result[4, :] = [self.indices[4], self.indices[5]]
        result[5, :] = [self.indices[5], self.indices[0]]
        return result


class Mesh:
    def __init__(self):
        self.elements = []
        self.ribs = []
        self.points = []

    def number_of_points(self):
        return len(self.points)

    def number_of_elements(self):
        return len(self.elements)

    def number_of_ribs(self):
        return len(self.ribs)

    def load_from_file(self, filename):
        mesh_file = open(filename, 'r')
        first_row_data = mesh_file.readline().split()
        nElements = int(first_row_data[0])
        nPoints = int(first_row_data[1])
        for i in range(nElements):
            indices_str = mesh_file.readline().split()
            indices = [int(index) for index in indices_str]
            self.elements.append(Element(indices))
        for i in range(nPoints):
            xy = mesh_file.readline().split()
            self.points.append(Point2D(np.float64(xy[0]), np.float64(xy[1])))
        self.ribs_init()

    def ribs_init(self):
        n = self.number_of_points()
        local_indices = [[i - 1, i] for i in range(6)]
        history = np.zeros(shape=[n, n], dtype=np.bool)
        for elem in self.elements:
            for loc_1, loc_2 in local_indices:
                glob_1, glob_2 = elem.indices[loc_1], elem.indices[loc_2]
                if not history[glob_1, glob_2]:
                    history[glob_1, glob_2] = True
                    history[glob_2, glob_1] = True
                    initial_length, _, _ = distance_between(self.points[glob_1], self.points[glob_2])
                    self.ribs.append(Rib([glob_1, glob_2], initial_length))

    def fi0(self, i, j, k):
        pi, pj, pk = self.points[i], self.points[j], self.points[k]
        fi0_ijk, _, _ = fi(pi.x, pj.x, pk.x, pi.y, pj.y, pk.y)
        return fi0_ijk

    def get_array_of_coordinates(self):
        coordinates = np.zeros(2*self.number_of_points())
        for i in range(self.number_of_points()):
            coordinates[2 * i] = self.points[i].x
            coordinates[2 * i + 1] = self.points[i].y
        return coordinates


def distance_between(p1, p2):
    distance = np.sqrt((p1.x - p2.x)*(p1.x - p2.x) +
                       (p1.y - p2.y) * (p1.y - p2.y))
    grad_x = [(p1.x - p2.x)/distance, (p2.x - p1.x)/distance]
    grad_y = [(p1.y - p2.y)/distance, (p2.y - p1.y)/distance]
    return distance, grad_x, grad_y


def cosfi_down(xi, xj, xk, yi, yj, yk):
    pi, pj, pk = Point2D(xi, yi), Point2D(xj, yj), Point2D(xk, yk)
    lij, lij_x, lij_y = distance_between(pi, pj)
    lik, lik_x, lik_y = distance_between(pi, pk)

    lij_xi, lij_xj = lij_x
    lij_yi, lij_yj = lij_y
    lik_xi, lik_xk = lik_x
    lik_yi, lik_yk = lik_y

    value = lij*lik
    grad_x = [lij_xi * lik + lij * lik_xi, lik*lij_xj, lij * lik_xk]
    grad_y = [lij_yi * lik + lij * lik_yi, lik * lij_yj, lij * lik_yk]
    return value, grad_x, grad_y


def cosfi_up(xi, xj, xk, yi, yj, yk):
    value = (xj - xi)*(xk - xi) + (yj - yi)*(yk - yi)
    grad_x = [2 * xi - xj - xk, xk - xi, xj - xi]
    grad_y = [2 * yi - yj - yk, yk - yi, yj - yi]
    return value, grad_x, grad_y


def cosfi(xi, xj, xk, yi, yj, yk):
    up, up_x, up_y = cosfi_up(xi, xj, xk, yi, yj, yk)
    down, down_x, down_y = cosfi_down(xi, xj, xk, yi, yj, yk)

    up_xi, up_xj, up_xk = up_x
    up_yi, up_yj, up_yk = up_y
    down_xi, down_xj, down_xk = down_x
    down_yi, down_yj, down_yk = down_y

    value = up/down
    grad_x = np.array([up_xi * down - up * down_xi,
              up_xj * down - up * down_xj,
              up_xk * down - up * down_xk]) / (down * down)
    grad_y = np.array([up_yi * down - up * down_yi,
              up_yj * down - up * down_yj,
              up_yk * down - up * down_yk]) / (down * down)

    return value, grad_x, grad_y


def fi(xi, xj, xk, yi, yj, yk):
    cos_fi, cos_fi_x, cos_fi_y = cosfi(xi, xj, xk, yi, yj, yk)
    # print(cos_fi)
    temp = -1.0/np.sqrt(1.0 - cos_fi*cos_fi)
    value = np.arccos(cos_fi)
    grad_x = temp * cos_fi_x
    grad_y = temp * cos_fi_y
    return value, grad_x, grad_y


def choose_coords(coordinates, index):
    global changeable, fixed_coordinates
    if changeable[index]:
        return coordinates[index]
    else:
        return fixed_coordinates[index]


def energy_angle(coords):
    value = 0.0
    grad = np.zeros(len(coords), dtype=np.float64)
    global mesh, changeable, fixed_coordinates
    for elem in mesh.elements:
        for angle in elem.angles():
            i, j, k = angle
            xi = choose_coords(coords, 2 * i)
            xj = choose_coords(coords, 2 * j)
            xk = choose_coords(coords, 2 * k)
            yi = choose_coords(coords, 2 * i + 1)
            yj = choose_coords(coords, 2 * j + 1)
            yk = choose_coords(coords, 2 * k + 1)
            fi_ijk, grad_x, grad_y = fi(xi, xj, xk, yi, yj, yk)
            temp = (fi_ijk - mesh.fi0(i, j, k))

            [grad[2 * i], grad[2 * j], grad[2 * k]] = \
                [grad[2 * i], grad[2 * j], grad[2 * k]] + 2.0 * temp * grad_x
            [grad[2 * i + 1], grad[2 * j + 1], grad[2 * k + 1]] = \
                [grad[2 * i + 1], grad[2 * j + 1], grad[2 * k + 1]] + 2.0 * temp * grad_y

            value += temp*temp
    return value, (grad*changeable)


def energy_length(coords):
    value = 0.0
    n = len(coords)
    grad = np.zeros(n, dtype=np.float64)
    global mesh, changeable, fixed_coordinates
    for rib in mesh.ribs:
        i, j = rib.indices
        xi = choose_coords(coords, 2 * i)
        xj = choose_coords(coords, 2 * j)
        yi = choose_coords(coords, 2 * i + 1)
        yj = choose_coords(coords, 2 * j + 1)
        l_ij, grad_x, grad_y = distance_between(Point2D(xi, yi), Point2D(xj, yj))
        l_ij_0 = rib.l0
        temp = l_ij - l_ij_0
        grad[2 * i] += 2.0 * temp * grad_x[0] / l_ij_0
        grad[2 * j] += 2.0 * temp * grad_x[1] / l_ij_0
        grad[2 * i + 1] += 2.0 * temp * grad_y[0] / l_ij_0
        grad[2 * j + 1] += 2.0 * temp * grad_y[1] / l_ij_0
        value += temp * temp / l_ij_0
    return value, (grad*changeable)


def full_energy(coords):
    w_angle = w_a
    w_length = w_l
    value_angle, grad_angle = energy_angle(coords)
    value_length, grad_length = energy_length(coords)
    return w_angle*value_angle + w_length*value_length, w_angle*grad_angle + w_length*grad_length


def get_constraints(filename, n_points):
    constr_file = open(filename, 'r')
    n_constraints = int(constr_file.readline())
    changeable = np.ones(shape=2*n_points, dtype=np.int8)
    fixed_coordinates = np.zeros(shape=2*n_points, dtype=np.float64)
    for i in range(n_constraints):
        data = constr_file.readline().split()
        index = int(data[0])
        mode = data[1]
        if mode == 'xy':
            changeable[2 * index] = 0
            changeable[2 * index + 1] = 0
            fixed_coordinates[2 * index] = np.float64(data[2])
            fixed_coordinates[2 * index + 1] = np.float64(data[3])
        if mode == 'x':
            changeable[2 * index] = 0
            fixed_coordinates[2 * index] = np.float64(data[2])
        if mode == 'y':
            changeable[2 * index + 1] = 0
            fixed_coordinates[2 * index + 1] = np.float64(data[2])
    return changeable, fixed_coordinates

def compute(mesh_filename, constr_filename, out_filename):
    mesh = Mesh()
    mesh.load_from_file(mesh_filename)
    changeable, fixed_coordinates = get_constraints(constr_filename, mesh.number_of_points())
    initial_coordinates = mesh.get_array_of_coordinates()
    coords = np.zeros(shape=2 * mesh.number_of_points(), dtype=np.float64)
    for i in range(mesh.number_of_points()):
        coords[2 * i] = choose_coords(initial_coordinates, 2 * i)
        coords[2 * i + 1] = choose_coords(initial_coordinates, 2 * i + 1)

    answer = opt.minimize(full_energy, jac=True, x0=coords, options={'maxiter': 10000, 'disp': True})
    out_file = open("result.txt", 'w')
    for i in range(mesh.number_of_points()):
        string = str(answer.x[2 * i]) + " " + str(answer.x[2 * i + 1])
        out_file.write(string + '\n')
    out_file.close()

    # analysis
    deformed = [Point2D(answer.x[2 * i], answer.x[2 * i + 1]) for i in range(mesh.number_of_points())]

    # distances before
    l_x_0 = (mesh.points[top[0]].x + mesh.points[top[1]].x) / 2. - (
                mesh.points[bottom[0]].x + mesh.points[bottom[1]].x) / 2.
    l_y_0 = mesh.points[right].y - mesh.points[left].y

    # distances after
    l_x = (deformed[top[0]].x + deformed[top[1]].x) / 2. - (deformed[bottom[0]].x + deformed[bottom[1]].x) / 2.
    l_y = deformed[right].y - deformed[left].y

    # deformations and poison ratio
    e_y = (l_y - l_y_0) / l_y_0
    e_x = (l_x - l_x_0) / l_x_0
    poison_ratio = - e_y / e_x
    # print("Poison ration is ", poison_ratio)

    return answer.fun, poison_ratio, e_x, e_y


if __name__ == "__main__":
    mesh_filename = 'mesh3x6.txt'
    constr_filename = 'constr3x6_21.txt'
    out_filename = "result.txt"
    print(compute(mesh_filename, constr_filename, out_filename))