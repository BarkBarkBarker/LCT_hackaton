import math
import pandas as pd
import numpy as np
from sympy import Plane, Point, Point3D, Line


class Tube:
    """
    Tube class, having information about its diam, connections with fittings
    """

    def __init__(self, type, id, diam, nominal_length, connection_shift):
        self.type = type
        self.id = id
        self.diam = diam
        self.nominal_length = nominal_length
        self.connection_shift = connection_shift
        self.n_nominal = None
        self.length = None
        self.in_port = None
        self.out_port = None

    def copy(self):
        return Tube(self.type, self.id, self.diam, self.nominal_length, self.connection_shift)


class Fitting:
    """
    Fitting class, having information about its diams, angles, number of ports, shifts of ports
    """

    def __init__(self, type, id, diams, angle, ports_shift):
        self.type = type
        self.id = id
        self.n_out_ports = len(diams) - 1
        self.diams = diams
        self.angle = angle
        self.ports_shift = ports_shift
        self.in_port = None
        self.out_ports = []
        self.reduction = None
        self.reduction_i = None

    def copy(self):
        return Fitting(self.type, self.id, self.diams, self.angle, self.ports_shift)


class FittingPort:
    """
    Port (hole) of fitting, must be connected with pipe
    """

    def __init__(self, point, diam):
        self.point = point
        self.diam = diam
        self.connected_tube = None
        self.owner = None
        self.reduction = False

    def copy(self):
        fit_copy = FittingPort(self.point, self.diam)
        fit_copy.connected_tube = None
        fit_copy.owner = None
        return fit_copy


class PlumbingUnit:
    """

    """

    def __init__(self, x, y, z, diameter):
        self.point = Point3D(x, y, z)
        self.diam = diameter
        self.wall = []
        self.wall_projection = None
        self.flat_point = None
        self.port = FittingPort(self.point, self.diam)

    def find_wall(self, walls):

        min_dist = 1e9
        for wall in walls:

            projection_point = wall.plane.projection(self.point)
            d = self.point.distance(projection_point)

            if d < min_dist and ((wall.l_point.x <= self.point.x <= wall.r_point.x) or (
                    wall.l_point.y >= self.point.y >= wall.r_point.y)):
                min_dist = d
                wall_result = wall
                wall_point = projection_point

        self.wall = wall_result
        self.wall_projection = wall_point
        wall_result.units.append(self)


class Wall:
    """
    self.l_point - most left or most upper point
    self.r_point - most right or most lower point
    """

    def __init__(self, x1, x2, y1, y2, width, height):
        if (x1 < x2 and y1 == y2) or (x1 == x2 and y1 > y2):
            self.l_point = Point(x1, y1, 0)
            self.r_point = Point(x2, y2, 0)
        elif (x1 > x2 and y1 == y2) or (x1 == x2 and y1 < y2):
            self.r_point = Point(x1, y1, 0)
            self.l_point = Point(x2, y2, 0)
        else:
            raise Exception('Incorrect wall coordinates')

        self.width = width
        self.height = height
        self.length = self.r_point.distance(self.l_point)

        self.plane = Plane(self.l_point, self.r_point, self.l_point + Point3D(0, 0, 1))
        self.base_line = Line(self.l_point, self.r_point)

        self.units = []

        self.delta_l = None

        self.neighbours = []

        self.l_neighbour = None
        self.r_neighbour = None

    pass


class Room:
    def __init__(self):
        self.walls = []
        self.fittings = []
        self.tubes = []
        self.fittings_allow = {}
        self.tubes_allow = {}
        self.water_raiser = None
        self.plumbing_units = []
        self.base_pipe_angle = 3  # degree
        self.continue_raiser = True  # raiser should be connected to vertical pipe up

    def add_wall(self, x1, x2, y1, y2, width, height):
        """
        Handler to Wall class to save in [walls]
        :return:
        """
        self.walls.append(Wall(x1, x2, y1, y2, width, height))

    def add_tube(self, port, point):
        """
        add tube between port and point, creates new port in point and return it
        :param port:
        :param point:
        :return:
        new_port
        """
        diam = port.diam
        length = point.distance(port.point)
        for tube in self.tubes_allow['Tube']:
            if tube.diam == diam:
                tube_copy = tube.copy()
                tube_copy.length = length + 2 * tube.connection_shift
                tube_copy.n_nominal = length / tube_copy.nominal_length
                tube_copy.in_port = port

                new_port = port.copy()
                new_port.point = point
                new_port.connected_tube = tube_copy
                new_port.wall = port.wall

                delta_point = new_port.point - port.point - (new_port.point.z - port.point.z) * Point3D(0, 0, 1)
                new_port.flat_point = port.flat_point + Point(delta_point.x + delta_point.y,
                                                              (new_port.point.z - port.point.z))

                port.connected_tube = tube_copy

                tube_copy.out_port = new_port

                self.tubes.append(tube_copy)
                return new_port

    def process_walls(self, units):
        """
        Finds connection between walls, delete unnecessary walls
        and calculate flat_points of units (their projections on walls)

        :return:
        array of corners flat coordinates (x)
        """
        # find connected walls
        for wall1 in self.walls:
            for wall2 in self.walls:
                if wall1 is not wall2:
                    if wall1.r_point == wall2.l_point or \
                            wall1.r_point == wall2.r_point or \
                            wall1.l_point == wall2.l_point:
                        if not (wall1 in wall2.neighbours) and not (wall2 in wall1.neighbours):
                            wall1.neighbours.append(wall2)
                            wall2.neighbours.append(wall1)

        edge_walls = []
        max_y = -1e9
        min_x = 1e9
        # find most top left edge wall
        for wall in self.walls:
            if len(wall.neighbours) == 1:
                if wall.l_point.x <= min_x:
                    if wall.l_point.y >= max_y:
                        max_y = wall.l_point.y
                        min_x = wall.l_point.x
                        top_left_edge_wall = wall

                edge_walls.append(wall)

        edge_walls.remove(top_left_edge_wall)

        # make a division into right and left neighbours
        cur_wall = top_left_edge_wall
        while not cur_wall == edge_walls[0]:
            cur_wall.r_neighbour = cur_wall.neighbours[0]
            cur_wall.neighbours[0].l_neighbour = cur_wall
            cur_wall.r_neighbour.neighbours.remove(cur_wall)
            cur_wall.neighbours.remove(cur_wall.r_neighbour)
            cur_wall = cur_wall.r_neighbour

        for plumb in units:
            plumb.find_wall(self.walls)

        # delete walls with no units connected on edge
        repeat_flag = True
        while repeat_flag:
            repeat_flag = False
            for wall in self.walls:
                if len(wall.units) == 0:
                    if wall.l_neighbour is None:
                        wall.r_neighbour.l_neighbour = None
                        self.walls.remove(wall)
                        repeat_flag = True
                    elif wall.r_neighbour is None:
                        wall.l_neighbour.r_neighbour = None
                        self.walls.remove(wall)
                        repeat_flag = True
        for wall in self.walls:
            if wall.l_neighbour is None:
                cur_wall = wall
                break
        delta_l = 0
        corners_x = []
        while cur_wall.r_neighbour is not None:  # move on neighbours to right edge wall
            cur_wall.delta_l = delta_l

            for unit in cur_wall.units:
                line = Line(cur_wall.l_point, cur_wall.l_point + Point3D(0, 0, cur_wall.height))

                l = line.distance(unit.wall_projection)
                unit.flat_point = Point(l + delta_l, unit.wall_projection.z)

            delta_l += cur_wall.length + cur_wall.r_neighbour.width / 2

            corners_x.append(delta_l)

            delta_l += cur_wall.width / 2

            cur_wall = cur_wall.r_neighbour
        else:  # most right wall
            cur_wall.delta_l = delta_l

            for unit in cur_wall.units:
                line = Line(cur_wall.l_point, cur_wall.l_point + Point3D(0, 0, cur_wall.height))

                l = line.distance(unit.wall_projection)

                unit.flat_point = Point(l + delta_l, unit.wall_projection.z)

        return corners_x

    def process_water_raiser(self, corners_x):
        """
        Separates plumbing_units into branches, add starting fitting
        corners_x - array of flat coordinates (x) of corners in room

        :return:
        port, branch array of units, branch corners for right and left sides
        """

        left_corners = []
        left_branch = []
        right_corners = []
        right_branch = []

        # separate branches and corners into branches
        for unit in self.plumbing_units:
            if unit.flat_point.x < self.water_raiser.flat_point.x:
                left_branch.append(unit)
            else:
                right_branch.append(unit)

        for corner in corners_x:
            if corner < self.water_raiser.flat_point.x:
                left_corners.append(corner)
            else:
                right_corners.append(corner)

        if len(left_branch) == 0 and not len(right_branch) == 0:  # have only right branch
            double = False
            branch = right_branch
            mirror = False
        elif len(right_branch) == 0 and not len(left_branch) == 0:  # have only left branch
            double = False
            branch = left_branch
            mirror = True
        elif not len(right_branch) == 0 and not len(left_branch) == 0:  # have both branches
            double = True
        else:
            raise Exception('No plumbing units found')

        if double:
            diams_branches = []
            angles_branches = []
            for branch in [left_branch, right_branch]:
                branch_diam = sorted(branch, key=lambda x: x.diam, reverse=True)

                diams_branches.append(branch_diam[0].diam)

                branch_x = sorted(branch, key=lambda x: x.flat_point.x)
                if self.water_raiser.flat_point.x > branch[0].flat_point.x:
                    branch_x = branch_x[::-1]

                angles = []
                for unit in branch_x:
                    angles.append(abs(math.atan(abs((self.water_raiser.flat_point.y - unit.flat_point.y) / (
                            self.water_raiser.flat_point.x - unit.flat_point.x)))))
                angles_branches.append(180 / math.pi * min(angles))

            if min(angles_branches) < self.base_pipe_angle:
                raise Exception(f'Cannot connect pipes with base angle={self.base_pipe_angle}, one of unit is too low')

            if self.continue_raiser:
                diams = [self.water_raiser.diam, self.water_raiser.diam] + diams_branches
                fitting = self.find_fitting('Quadruple', diams, 90 - self.base_pipe_angle)
                if diams_branches[0] < diams_branches[1]:
                    out_ports = self.install_fitting(fitting, self.water_raiser.port, self.water_raiser.wall, 180, 180,
                                                     0, '')
                else:
                    out_ports = self.install_fitting(fitting, self.water_raiser.port, self.water_raiser.wall, 0, 0, 0,
                                                     '')
                right_port = out_ports[1]
                left_port = out_ports[2]
            else:
                diams = [self.water_raiser.diam] + diams_branches
                fitting = self.find_fitting('Triple', diams, 90 - self.base_pipe_angle)
                if diams_branches[0] < diams_branches[1]:
                    out_ports = self.install_fitting(fitting, self.water_raiser.port, self.water_raiser.wall, 180, 180,
                                                     0, '')
                else:
                    out_ports = self.install_fitting(fitting, self.water_raiser.port, self.water_raiser.wall, 0, 0, 0,
                                                     '')
                right_port = out_ports[0]
                left_port = out_ports[1]
        else:
            branch_diam = sorted(branch, key=lambda x: x.diam, reverse=True)

            max_diam_branch = branch_diam[0].diam

            branch_x = sorted(branch, key=lambda x: x.flat_point.x)
            if self.water_raiser.flat_point.x > branch[0].flat_point.x:
                branch_x = branch_x[::-1]

            angles = []
            for unit in branch_x:
                angles.append(abs(math.atan(abs((self.water_raiser.flat_point.y - unit.flat_point.y) / (
                        self.water_raiser.flat_point.x - unit.flat_point.x)))))
            min_angle = 180 / math.pi * min(angles)
            if min_angle < self.base_pipe_angle:
                raise Exception(f'Cannot connect pipes with base angle={self.base_pipe_angle}, one of unit is too low')

            if self.continue_raiser:
                diams = [self.water_raiser.diam, self.water_raiser.diam, max_diam_branch]
                fitting = self.find_fitting('Triple', diams, 90 - self.base_pipe_angle)
                if mirror:
                    out_ports = self.install_fitting(fitting, self.water_raiser.port, self.water_raiser.wall, 180, 180,
                                                     0, '')
                else:
                    out_ports = self.install_fitting(fitting, self.water_raiser.port, self.water_raiser.wall, 0, 0, 0,
                                                     '')

                out_port = out_ports[1]
            else:
                diams = [self.water_raiser.diam, max_diam_branch]
                fitting = self.find_fitting('Turning', diams, 90 - self.base_pipe_angle)
                if mirror:
                    out_ports = self.install_fitting(fitting, self.water_raiser.port, self.water_raiser.wall, 180, 180,
                                                     0, '')
                else:
                    out_ports = self.install_fitting(fitting, self.water_raiser.port, self.water_raiser.wall, 0, 0, 0,
                                                     '')

                out_port = out_ports[0]

            if len(left_branch) == 0 and not len(right_branch) == 0:
                left_port = None
                right_port = out_port

                line = Line(self.water_raiser.wall.l_point,
                            self.water_raiser.wall.l_point + Point3D(0, 0, self.water_raiser.wall.height))
                projection = self.water_raiser.wall.plane.projection(right_port.point)
                l = line.distance(projection)
                right_port.flat_point = Point(l + self.water_raiser.wall.delta_l, projection.z)
            else:
                left_port = out_port
                right_port = None

                line = Line(self.water_raiser.wall.l_point,
                            self.water_raiser.wall.l_point + Point3D(0, 0, self.water_raiser.wall.height))
                projection = self.water_raiser.wall.plane.projection(left_port.point)
                l = line.distance(projection)
                left_port.flat_point = Point(l + self.water_raiser.wall.delta_l, projection.z)

        return left_port, left_branch, left_corners, right_port, right_branch, right_corners

    def process_plumbing_branch(self, in_port, units, corners, direction):
        """
        Process connections of plumbing units on one direction from water_raiser
            in_port - port of water_raiser to this branch
            units - units in this branch
            corners - corners flat coordinates in this branch
            direction - 'left' or 'right'
        :return:
        """

        if len(units) == 0:
            return

        branch_x = sorted(units, key=lambda x: x.flat_point.x)
        if direction == 'left':
            branch_x = branch_x[::-1]

        cur_port = in_port

        for i in range(len(branch_x)):
            for corner_x in corners:  # case when between two ports had corner
                if (cur_port.flat_point.x < corner_x < branch_x[i].flat_point.x) or (
                        cur_port.flat_point.x > corner_x > branch_x[i].flat_point.x):
                    corner_y = cur_port.flat_point.y + (math.tan(self.base_pipe_angle / 180 * math.pi)) * abs(
                        cur_port.flat_point.x - corner_x)
                    corner_point = cur_port.point + corner_y * Point3D(0, 0, 1)
                    corner_point += (
                                            cur_port.owner.wall.r_point - cur_port.owner.wall.l_point) / cur_port.owner.wall.length * (
                                            cur_port.flat_point.x - corner_x)

                    new_port = self.add_tube(cur_port, corner_point)

                    turn45 = self.find_fitting('Turning', [cur_port.diam, cur_port.diam], 45)
                    ports = self.install_fitting(turn45, new_port, cur_port.wall, (90 - self.base_pipe_angle),
                                                 -90, 0, direction + ' corner 1')
                    new_port.owner = ports[0].owner
                    new_port = ports[0]
                    turn45 = turn45.copy()
                    ports = self.install_fitting(turn45, new_port, cur_port.wall, (90 - self.base_pipe_angle),
                                                 -90, -45, direction + ' corner 2')
                    cur_port = ports[0]

            branch_y = cur_port.flat_point.y + math.tan(self.base_pipe_angle / 180 * math.pi) * abs(
                cur_port.flat_point.x - branch_x[i].flat_point.x)
            if not (i == len(branch_x) - 1):
                point = (branch_x[i].wall.r_point - branch_x[i].wall.l_point) / branch_x[i].wall.length * (
                        branch_x[i].flat_point.x - cur_port.flat_point.x)
                point += Point(0, 0, 1) * (branch_y - cur_port.flat_point.y)

                port = self.add_tube(cur_port, point)

                outer_branch_diams = [branch_x[j].port.diam for j in range(i + 1, len(branch_x))]

                diams = [cur_port.diam, max(outer_branch_diams), branch_x[i].port.diam]

                fitting = self.find_fitting('Triple', diams, 45)

                ports = self.install_fitting(fitting, port, cur_port.wall, (90 - self.base_pipe_angle), 0, 0, '')

                port.owner = ports[0].owner

                new_port = ports[1]
                cur_port = ports[0]

                fitting = self.find_fitting('Turning', [new_port.diam, new_port.diam], 45)

                ports = self.install_fitting(fitting, new_port, cur_port.wall, (45 - self.base_pipe_angle), 0, 0, '')

                port = self.add_tube(ports[0], branch_x[i].port.point)

                fitting = self.find_fitting('Turning', [new_port.diam, new_port.diam], 87)

                ports = self.install_fitting(fitting, port, cur_port.wall, (-self.base_pipe_angle), 0, -90, '')

                port.owner = ports[0].owner

                branch_x[i].port = ports[0]
            else:
                point = (branch_x[i].wall.r_point - branch_x[i].wall.l_point) / branch_x[i].wall.length * (
                        branch_x[i].flat_point.x - cur_port.flat_point.x)
                point += Point(0, 0, 1) * (branch_y - cur_port.flat_point.y)

                diams = [cur_port.diam, branch_x[i].port.diam]

                fitting = self.find_fitting('Turning', diams, 45)

                ports = self.install_fitting(fitting, cur_port, cur_port.wall, (90 - self.base_pipe_angle), 0, 0, '')

                new_port = ports[0]

                fitting = self.find_fitting('Turning', [new_port.diam, new_port.diam], 45)

                ports = self.install_fitting(fitting, new_port, cur_port.wall, (45 - self.base_pipe_angle), 0, 0, '')

                port = self.add_tube(ports[0], branch_x[i].port.point)

                fitting = self.find_fitting('Turning', [new_port.diam, new_port.diam], 87)

                ports = self.install_fitting(fitting, port, cur_port.wall, (-self.base_pipe_angle), 0, -90, '')

                port.owner = ports[0].owner

                branch_x[i].port = ports[0]

    def find_fitting(self, type, diams, angle):
        """
        Finds in self.fittings_allow suitable fittings, returns this fitting_model
        :return:
        fitting_model or None if there's no such fitting
        """
        for fitting in self.fittings_allow[type]:
            if fitting.diams == diams and fitting.angle == angle:
                result = fitting
                break
        else:
            for fitting in self.fittings_allow[type]:
                if fitting.angle == angle:
                    for i in range(1, len(diams)):
                        if not diams[i] == fitting.diams[i]:
                            for reduction in self.fittings_allow['Reduction']:
                                if reduction.diams[0] == fitting.diams[i] and reduction.diams[1] == diams[i]:
                                    if not i == len(diams) - 1:
                                        if fitting.diams[0:i] + [reduction.diams[1]] + fitting.diams[i + 1:] == diams:
                                            result = fitting.copy()
                                            result.reduction_i = i
                                            result.reduction = reduction
                                            return result
                                    else:
                                        if fitting.diams[0:i] + [reduction.diams[1]] == diams:
                                            result = fitting.copy()
                                            result.reduction_i = i
                                            result.reduction = reduction
                                            return result

            raise Exception(
                f"Can't find fitting of type={type} with diameters {diams} and angle = {angle}°")
        return result

    def install_fitting(self, fitting_model, in_port, wall, rotation_local_x, rotation_local_y, rotation_local_z, type):
        """
        Creates Fitting in plane of wall
        with rotation_plane (polar angle in radians) and rotation_deep (angle to normal of wall)
        :return:
        array of out_ports
        """
        angle1 = rotation_local_x / 180 * math.pi
        angle2 = rotation_local_y / 180 * math.pi
        angle3 = rotation_local_z / 180 * math.pi

        if fitting_model.reduction is not None:
            fitting_model.ports_shift[fitting_model.reduction_i] += fitting_model.reduction.ports_shift[0]

        rotx = lambda alpha: np.array([[1, 0, 0],
                                       [0, np.cos(alpha), -np.sin(alpha)],
                                       [0, np.sin(alpha), np.cos(alpha)]])

        roty = lambda alpha: np.array([[np.cos(alpha), 0, np.sin(alpha)],
                                       [0, 1, 0],
                                       [-np.sin(alpha), 0, np.cos(alpha)]])

        rotz = lambda alpha: np.array([[np.cos(alpha), -np.sin(alpha), 0],
                                       [np.sin(alpha), np.cos(alpha), 0],
                                       [0, 0, 1]])

        if Line.is_parallel(wall.base_line, Line((0, 0, 0), (1, 0, 0))):  # wall lay on x-axis
            matrix_rot11 = roty(angle1)  # rotate around y-axis
            matrix_rot12 = roty(math.pi / 2 - angle1)  # rotate around y-axis
            matrix_rot2 = rotx(angle2)  # rotate around x-axis
        elif Line.is_parallel(wall.base_line, Line((0, 0, 0), (0, 1, 0))):  # wall lay on y-axis
            matrix_rot11 = rotx(angle1)  # rotate around x-axis
            matrix_rot12 = rotx(math.pi / 2 - angle1)  # rotate around x-axis
            matrix_rot2 = roty(angle2)  # rotate around y-axis
        else:
            raise Exception('Walls must be parallel to x- or y-axis')

        matrix_rot3 = rotz(angle3)  # rotate around z-axis

        fitting_copy = fitting_model.copy()

        if type == 'right corner 1':
            matrix = np.matmul(matrix_rot12, np.matmul(rotz(math.pi / 2), matrix_rot2))
        elif type == 'right corner 2':
            matrix = np.matmul(rotz(-math.pi / 4), np.matmul(matrix_rot12, np.matmul(rotz(math.pi / 2), matrix_rot2)))
        elif type == 'left corner 1':
            matrix = np.matmul(matrix_rot12, np.matmul(rotz(-math.pi / 2), np.matmul(matrix_rot2, rotz(math.pi))))
        elif type == 'left corner 2':
            matrix = np.matmul(rotz(math.pi / 4), np.matmul(matrix_rot12, np.matmul(rotz(-math.pi / 2),
                                                                                    np.matmul(matrix_rot2,
                                                                                              rotz(math.pi / 2)))))
        else:
            matrix = np.matmul(matrix_rot11, np.matmul(matrix_rot2, matrix_rot3))

        for i in range(fitting_copy.n_out_ports):
            point_delta = fitting_copy.ports_shift[i]
            point_vector = np.array(point_delta)
            rotated_vector = np.matmul(matrix, point_vector)
            fitting_copy.ports_shift[i] = Point3D(rotated_vector)

        fitting_copy.in_port = in_port

        line = Line(wall.l_point, wall.l_point + Point3D(0, 0, wall.height))

        for i in range(fitting_model.n_out_ports):
            port_copy = in_port.copy()
            port_copy.owner = fitting_copy
            port_copy.diam = fitting_copy.diams[1 + i]
            port_copy.point += fitting_copy.ports_shift[i]

            line = Line(wall.l_point, wall.l_point + Point3D(0, 0, wall.height))
            projection = wall.plane.projection(port_copy.point)
            l = line.distance(projection)
            port_copy.flat_point = Point(l + wall.delta_l, projection.z)

            fitting_copy.out_ports.append(port_copy)

            fitting_copy.wall = wall
            port_copy.wall = wall

        if fitting_model.reduction is not None:
            fitting_copy.out_ports[fitting_model.reduction_i - 1].reduction = True
            fitting_copy.out_ports[fitting_model.reduction_i - 1].diam = fitting_model.reduction.diams[1]
            fitting_copy.diams[fitting_model.reduction_i] = fitting_model.reduction.diams[1]
            fitting_copy.reduction = fitting_model.reduction

        self.fittings.append(fitting_copy)
        return fitting_copy.out_ports

    def export_data(self):
        """
        Exports information about pipes, fittings, it's connection
        :return:
        """
        type = []
        id_table = []
        diams = []
        x = []
        y = []
        z = []
        x1 = []
        x2 = []
        y1 = []
        y2 = []
        z1 = []
        z2 = []
        connections = []
        connections_diam = []
        lengths = []

        units = [self.water_raiser] + self.plumbing_units

        for i in range(len(units)):
            type.append('Сантех.узел')
            id_table.append(i)
            diams.append(units[i].port.diam)
            x.append(round(units[i].port.point.x, 3))
            y.append(round(units[i].port.point.y, 3))
            z.append(round(units[i].port.point.z, 3))
            N = i
            units[i].id_table = id_table[-1]
        for j in range(len(self.fittings)):
            type.append(f'Фиттинг {self.fittings[j].type} с ИД {self.fittings[j].id}')
            id_table.append(N + j + 1)
            diams.append(' '.join([str(elem) for elem in self.fittings[j].diams]))
            x.append(round(self.fittings[j].in_port.point.x, 3))
            y.append(round(self.fittings[j].in_port.point.y, 3))
            z.append(round(self.fittings[j].in_port.point.z, 3))
            self.fittings[j].id_table = id_table[-1]
            if self.fittings[j].reduction is not None:
                type.append(self.fittings[j].reduction.type)
                N += 1
                id_table.append(N + j + 1)
                diams.append(self.fittings[j].reduction.diams)
                x.append(x[-1])
                y.append(y[-1])
                z.append(z[-1])
                x1.append(' ')
                x2.append(' ')
                y1.append(' ')
                y2.append(' ')
                z1.append(' ')
                z2.append(' ')

                text = str(id_table[-2]) + ' - ' + str(id_table[-1])
                connections.append(text)
                connections_diam.append(self.fittings[j].reduction.diams)
                lengths.append(' ')

        for tube in self.tubes:
            text = str(tube.in_port.owner.id_table) + ' - ' + str(tube.out_port.owner.id_table)
            connections.append(text)
            connections_diam.append(tube.diam)
            lengths.append(round(tube.length, 3))
            x1.append(round(tube.in_port.point.x, 3))
            x2.append(round(tube.out_port.point.x, 3))
            y1.append(round(tube.in_port.point.y, 3))
            y2.append(round(tube.out_port.point.y, 3))
            z1.append(round(tube.in_port.point.z, 3))
            z2.append(round(tube.out_port.point.z, 3))

        data1 = {'ИД': id_table, 'Тип': type, 'Диаметр': diams, 'X': x, 'Y': y, 'Z': z}
        data2 = {'Звено': connections, 'Диаметр': connections_diam, 'Длина [м]': lengths, 'X1': x1, 'X2': x2, 'Y1': y1,
                 'Y2': y2, 'Z1': z1, 'Z2': z2}

        df1 = pd.DataFrame(data1)
        df2 = pd.DataFrame(data2)

        with pd.ExcelWriter('dataframe.xlsx') as writer:
            df1.to_excel(writer, sheet_name='Узлы')
            df2.to_excel(writer, sheet_name='Трубы')

    def calculate_tubes(self):
        """
        Main function to run calculations
        :return:
        """
        units = [self.water_raiser] + self.plumbing_units

        corners = self.process_walls(units)

        l_port, l_branch, l_corners, r_port, r_branch, r_corners = self.process_water_raiser(corners)

        self.process_plumbing_branch(l_port, l_branch, l_corners, 'left')
        self.process_plumbing_branch(r_port, r_branch, r_corners, 'right')

        self.export_data()


room = Room()

fittings_allow = {'Triple': [Fitting('Triple', 101, [110, 110, 110], 87,
                                     [Point3D(0, 0, 0.015), Point3D(0.005, 0, 0.005)]),
                             Fitting('Triple', 102, [110, 50, 110], 87,
                                     [Point3D(0, 0, 0.015), Point3D(0.005, 0, 0.005)]),
                             Fitting('Triple', 103, [50, 50, 50], 87,
                                     [Point3D(0, 0, 0.015), Point3D(0.005, 0, 0.005)]),
                             Fitting('Triple', 104, [110, 110, 110], 45,
                                     [Point3D(0, 0, 0.015), Point3D(0.005, 0, 0.005)]),
                             Fitting('Triple', 106, [50, 50, 50], 45,
                                     [Point3D(0, 0, 0.015), Point3D(0.005, 0, 0.005)])
                             ],

                  'Quadruple': [Fitting('Quadruple', 201, [110, 110, 110, 50], 87,
                                        [Point3D(0, 0, 0.015), Point3D(0.005, 0, 0.005), Point3D(-0.005, 0, 0.005)])],

                  'Turning': [Fitting('Turning', 310, [40, 40], 87, [Point3D(0.005, 0, 0.005)]),
                              Fitting('Turning', 305, [50, 50], 87, [Point3D(0.005, 0, 0.005)]),
                              Fitting('Turning', 306, [110, 110], 87, [Point3D(0.005, 0, 0.005)]),
                              Fitting('Turning', 309, [40, 40], 45, [Point3D(0.005, 0, 0.005)]),
                              Fitting('Turning', 304, [50, 50], 45, [Point3D(0.005, 0, 0.005)]),
                              Fitting('Turning', 300, [110, 110], 45, [Point3D(0.005, 0, 0.005)])],

                  'Reduction': [Fitting('Reduction', 401, [110, 50], 0, [Point3D(0.002, 0, 0.015)]),
                                Fitting('Reduction', 402, [50, 40], 0, [Point3D(0.002, 0, 0.015)])]}

tubes_allow = {'Tube': [Tube('Tube', 500, 110, 3, 0.06),
                        Tube('Tube', 501, 50, 3, 0.044),
                        Tube('Tube', 502, 40, 2, 0.039)],
               'Nozzle': []}

room.fittings_allow = fittings_allow
room.tubes_allow = tubes_allow

room.add_wall(-0.15, 1.568, -0.1535, -0.1535, 0.2535, 1.2)
room.add_wall(1.568, 1.568, -0.1535, -2.310, 0.1, 2)
room.add_wall(1.568, -0.15, -2.310, -2.310, 0.05, 2)
room.add_wall(-0.15, -0.15, -0.1535, -1, 0.1, 2)

raiser = PlumbingUnit(0, 0, 0, 110)
toilet = PlumbingUnit(0.34745, -0.153, 0.2485, 110)
sink = PlumbingUnit(1.568, -1.2301, 0.616, 50)
bath = PlumbingUnit(1.568, -1.943, 0.19850, 50)
boss = PlumbingUnit(1.585069, -2.648512, 0.410237, 40)

room.water_raiser = raiser
room.plumbing_units = [toilet, sink, bath]

room.calculate_tubes()

room.export_data()

print(1)
