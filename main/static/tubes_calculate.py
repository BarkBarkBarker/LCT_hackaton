import math
import pandas as pd
import numpy as np
from sympy import Plane, Point, Point3D, Line
import sqlite3


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
        shift_copy = []
        for i in self.ports_shift:
            shift_copy.append(i)
        return Fitting(self.type, self.id, self.diams, self.angle, shift_copy)

    def a(self):
        for port in self.out_ports:
            print(port.point.x.evalf(), ',', port.point.y.evalf(), ',', port.point.z.evalf())
        print(self.in_port.point.x.evalf(), ',',self.in_port.point.y.evalf(), ',',self.in_port.point.z.evalf())

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
        self.wall = None
        self.wall_projection = None
        self.flat_point = None
        self.room = None
        self.direct = False
        self.port = FittingPort(self.point, self.diam)

    def find_wall(self, walls):

        wall_result = None
        wall_point = None

        min_dist = 1e9
        for wall in walls:

            projection_point = wall.plane.projection(self.point)
            d = self.point.distance(projection_point)
            if d == 0:
                raise(Exception('Plumbing unit point lays on wall plane'))
            if d < min_dist and ((wall.l_point.x <= self.point.x <= wall.r_point.x) or (
                    wall.l_point.y >= self.point.y >= wall.r_point.y)):
                min_dist = d
                wall_result = wall
                wall_point = projection_point

        if wall_result is None and self is self.room.water_raiser:
            for unit in self.room.plumbing_units:
                if unit.diam == 110:
                    wall_point = unit.wall.plane.projection(self.point)
                    wall_result = unit.wall
                    break

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
        self.double = False

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
                            wall1.l_point == wall2.l_point or \
                            wall1.l_point == wall2.r_point:
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

            # changed = False
            # for p1 in [cur_wall.l_point, cur_wall.r_point]:
            #     for p2 in [cur_wall.r_neighbour.l_point, cur_wall.r_neighbour.r_point]:
            #         if p1 == p2 and not changed:
            #             arr = [cur_wall.l_point, cur_wall.r_point]
            #             cur_wall.r_point = p1
            #             arr.remove(p1)
            #             cur_wall.l_point = arr[0]
            #
            #             arr = [cur_wall.r_neighbour.l_point, cur_wall.r_neighbour.r_point]
            #             cur_wall.r_neighbour.l_point = p2
            #             arr.remove(p2)
            #             cur_wall.r_neighbour.r_point = arr[0]
            #             changed = True

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

        line_direction = Line(self.water_raiser.wall.r_point, self.water_raiser.wall.l_point)

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
            mirror = True
        elif len(right_branch) == 0 and not len(left_branch) == 0:  # have only left branch
            double = False
            branch = left_branch
            mirror = False
        elif not len(right_branch) == 0 and not len(left_branch) == 0:  # have both branches
            double = True
        else:
            raise Exception('No plumbing units found')

        if double:  # add quadruple fitting
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

            diams = [self.water_raiser.diam, self.water_raiser.diam] + diams_branches
            fitting = self.find_fitting('Quadruple', diams, 90 - self.base_pipe_angle)
            if diams_branches[0] < diams_branches[1]:  # left one has lower diameter then rotate 180deg on Oz
                out_ports = self.install_fitting(fitting, self.water_raiser.port, self.water_raiser.wall, line_direction, 0, 0, 0, '')
            else:
                out_ports = self.install_fitting(fitting, self.water_raiser.port, self.water_raiser.wall, line_direction, 0, 0, 180, '')
            right_port = out_ports[1]
            left_port = out_ports[2]

        else:  # add triple fitting
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

            diams = [self.water_raiser.diam, self.water_raiser.diam, max_diam_branch]
            fitting = self.find_fitting('Triple', diams, 90 - self.base_pipe_angle)
            if mirror:
                out_ports = self.install_fitting(fitting, self.water_raiser.port, self.water_raiser.wall, line_direction, 0, 0, 0, '')
            else:
                out_ports = self.install_fitting(fitting, self.water_raiser.port, self.water_raiser.wall, line_direction, 0, 0, 180, '')

            out_port = out_ports[1]

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
        line_direction = Line(self.water_raiser.wall.l_point, self.water_raiser.wall.r_point)

        cur_port = in_port
        cur_wall = self.water_raiser.wall

        for i in range(len(branch_x)):
            for corner_x in corners:  # case when between two ports had corner
                if (cur_port.flat_point.x < corner_x < branch_x[i].flat_point.x) or (
                        cur_port.flat_point.x > corner_x > branch_x[i].flat_point.x):

                    if direction == 'left':
                        next_wall = cur_port.wall.l_neighbour
                    else:
                        next_wall = cur_port.wall.r_neighbour

                    corner_y = (math.tan(self.base_pipe_angle / 180 * math.pi)) * abs(
                        cur_port.flat_point.x - corner_x)
                    corner_point = cur_port.point + corner_y * Point3D(0, 0, 1)

                    turn45 = self.find_fitting('Turning', [cur_port.diam, cur_port.diam], 45)

                    delta_tube_length = (turn45.ports_shift[0].z + turn45.ports_shift[0].x)

                    corner_point += (cur_port.owner.wall.r_point - cur_port.owner.wall.l_point) / cur_port.owner.wall.length * (
                                            corner_x - cur_port.flat_point.x - (corner_x - cur_port.flat_point.x)/abs(corner_x - cur_port.flat_point.x)*delta_tube_length*np.cos(self.base_pipe_angle/180*np.pi))

                    corner_point += (cur_port.owner.wall.r_point - cur_port.owner.wall.l_point) / cur_port.owner.wall.length * next_wall.width

                    new_port = self.add_tube(cur_port, corner_point)

                    ports = self.install_fitting(turn45, new_port, cur_wall, line_direction, 0, -(90 - self.base_pipe_angle), 90, direction)
                    new_port.owner = ports[0].owner
                    new_port = ports[0]

                    if direction == 'left':
                        new_port.flat_point += delta_tube_length * Point(-np.cos(self.base_pipe_angle / 180 * np.pi),
                                                                         np.sin(self.base_pipe_angle / 180 * np.pi))
                    else:
                        new_port.flat_point += delta_tube_length * Point(np.cos(self.base_pipe_angle / 180 * np.pi),
                                                                         np.sin(self.base_pipe_angle / 180 * np.pi))

                    turn45 = turn45.copy()
                    ports = self.install_fitting(turn45, new_port, cur_wall, line_direction, -(90 - self.base_pipe_angle), -45, 0, direction)
                    if direction == 'left':
                        ports[0].flat_point += delta_tube_length * Point(-np.cos(self.base_pipe_angle / 180 * np.pi),
                                                                         np.sin(self.base_pipe_angle / 180 * np.pi))
                    else:
                        ports[0].flat_point += delta_tube_length * Point(np.cos(self.base_pipe_angle / 180 * np.pi),
                                                                         np.sin(self.base_pipe_angle / 180 * np.pi))
                    if direction == 'left':
                        for p in [cur_port.wall.l_neighbour.l_point, cur_port.wall.l_neighbour.r_point]:
                            if line_direction.p1 != p:
                                line_direction = Line(p, line_direction.p1)
                                break
                    else:
                        for p in [cur_port.wall.r_neighbour.l_point, cur_port.wall.r_neighbour.r_point]:
                            if line_direction.p2 != p:
                                line_direction = Line(line_direction.p2, p)
                                break
                    cur_wall = next_wall
                    cur_port = ports[0]

            branch_y = cur_port.flat_point.y + math.tan(self.base_pipe_angle / 180 * math.pi) * abs(
                cur_port.flat_point.x - branch_x[i].flat_point.x)

            if not (i == len(branch_x) - 1):
                point = cur_port.point
                if direction == 'right':
                    point += (branch_x[i].wall.r_point - branch_x[i].wall.l_point) / branch_x[i].wall.length * abs(
                            branch_x[i].flat_point.x - cur_port.flat_point.x)  # delta x or y
                else:
                    point += (branch_x[i].wall.l_point - branch_x[i].wall.r_point) / branch_x[i].wall.length * abs(
                        branch_x[i].flat_point.x - cur_port.flat_point.x)  # delta x or y
                point += Point(0, 0, 1) * (branch_y - cur_port.flat_point.y)  # delta z

                outer_branch_diams = [branch_x[j].port.diam for j in range(i + 1, len(branch_x))]

                diams = [cur_port.diam, max(outer_branch_diams), branch_x[i].port.diam]

                fitting_short_path = self.find_fitting('Triple', diams, 87)
                fitting_short_path_turn = self.find_fitting('Turning', [branch_x[i].port.diam, branch_x[i].port.diam], 45)
                height_lim = fitting_short_path.ports_shift[1].x

                #dp = (branch_x[i].wall_projection - branch_x[i].point)/branch_x[i].point.distance(branch_x[i].wall_projection)
                #rotated_port_point = point + dp * height_lim + Point3D(0, 0, 1) * height_lim

                if point.distance(branch_x[i].port.point) - height_lim <= fitting_short_path_turn.ports_shift[0].distance(Point3D(0, 0, 0)):
                    point -= (point - cur_port.point) / point.distance(cur_port.point) * fitting_short_path.ports_shift[1].z
                    port = self.add_tube(cur_port, point)
                    ports = self.install_fitting(fitting_short_path, port, cur_wall, line_direction, 0, -(90 - self.base_pipe_angle), 45, '')
                    port.owner = ports[0].owner

                    new_port = ports[1]
                    cur_port = ports[0]

                    ports = self.install_fitting(fitting_short_path_turn, new_port, cur_wall, line_direction, -45, 0, 90, '')

                    branch_x[i].port = ports[0]
                else:
                    fitting_long_path = self.find_fitting('Triple', diams, 45)

                    point -= (point - cur_port.point)/point.distance(cur_port.point) * fitting_long_path.ports_shift[1].z

                    port = self.add_tube(cur_port, point)

                    ports = self.install_fitting(fitting_long_path, port, cur_wall, line_direction, 0, -(90 - self.base_pipe_angle), 0, '')

                    port.owner = ports[0].owner

                    new_port = ports[1]
                    cur_port = ports[0]

                    fitting1 = self.find_fitting('Turning', [new_port.diam, new_port.diam], 45)
                    fitting2 = self.find_fitting('Turning', [new_port.diam, new_port.diam], 87)

                    ports = self.install_fitting(fitting1, new_port, cur_wall, line_direction, 0, -(45 - self.base_pipe_angle), 0, '')

                    point = branch_x[i].point - (branch_x[i].point - ports[0].point) / branch_x[i].point.distance(ports[0].point) * fitting2.ports_shift[0].z

                    port = self.add_tube(ports[0], point)

                    ports = self.install_fitting(fitting2, port, cur_wall, line_direction, 0, (self.base_pipe_angle), 90, '')

                    port.owner = ports[0].owner

                    branch_x[i].port = ports[0]
            else:
                if not branch_x[i].direct:

                    diams = [cur_port.diam, branch_x[i].port.diam]

                    fitting45 = self.find_fitting('Turning', diams, 45)

                    point = cur_port.point

                    point += (branch_x[i].wall.r_point - branch_x[i].wall.l_point) / branch_x[i].wall.length * abs(branch_x[i].flat_point.x - cur_port.flat_point.x - fitting45.ports_shift[0].z * np.cos(self.base_pipe_angle/180*np.pi))  # delta x or y
                    point += Point(0, 0, 1) * (branch_y - cur_port.flat_point.y - fitting45.ports_shift[0].z * np.sin(self.base_pipe_angle/180*np.pi))  # delta z

                    port = self.add_tube(cur_port, point)

                    ports = self.install_fitting(fitting45, port, cur_wall, line_direction, 0, -(90 - self.base_pipe_angle), 0, '')

                    port.owner = ports[0].owner

                    new_port = ports[0]

                    fitting45 = fitting45.copy()

                    ports = self.install_fitting(fitting45, new_port, cur_wall, line_direction, 0, -(45 - self.base_pipe_angle), 0, '')

                    fitting87 = self.find_fitting('Turning', diams, 87)

                    point = branch_x[i].point - (branch_x[i].point - ports[0].point) / branch_x[i].point.distance(
                        ports[0].point) * fitting87.ports_shift[0].z

                    port = self.add_tube(ports[0], point)

                    ports = self.install_fitting(fitting87, port, cur_wall, line_direction, 0, (self.base_pipe_angle), 90, '')

                    port.owner = ports[0].owner

                    branch_x[i].port = ports[0]
                else:
                    port = self.add_tube(cur_port, branch_x[i].port.point)
                    branch_x[i].port = port
                    branch_x[i].port.owner = branch_x[i]

    def find_fitting(self, type, diams, angle):
        """
        Finds in self.fittings_allow suitable fittings, returns this fitting_model
        :return:
        fitting_model or None if there's no such fitting
        """
        for fitting in self.fittings_allow[type]:
            if fitting.diams == diams and fitting.angle == angle:
                result = fitting.copy()
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
            if type == 'Quadruple' and diams[2] != diams[3]:
                diams[2], diams[3] = diams[3], diams[2]
                for fitting in self.fittings_allow[type]:
                    if fitting.angle == angle:
                        if fitting.diams == diams:
                            result = fitting.copy()
                            return result
                        for i in range(1, len(diams)):
                            if not diams[i] == fitting.diams[i]:
                                for reduction in self.fittings_allow['Reduction']:
                                    if reduction.diams[0] == fitting.diams[i] and reduction.diams[1] == diams[i]:
                                        if not i == len(diams) - 1:
                                            if fitting.diams[0:i] + [reduction.diams[1]] + fitting.diams[
                                                                                           i + 1:] == diams:
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

    def install_fitting(self, fitting_model, in_port, wall, line_direction, rotation_local_x, rotation_local_y, rotation_local_z, direction):
        """
        Creates Fitting in plane of wall
        with rotation_plane (polar angle in radians) and rotation_deep (angle to normal of wall)
        :return:
        array of out_ports
        """

        base_angle = float(Line(Point3D(0, 0, 0), Point3D(-1, 0, 0)).angle_between(line_direction).evalf())

        angle1 = rotation_local_x / 180 * math.pi
        angle2 = rotation_local_y / 180 * math.pi
        angle3 = rotation_local_z / 180 * math.pi

        if fitting_model.reduction is not None and fitting_model.reduction_i == 0:
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

        base_mat = rotz(base_angle)
        matx = rotx(angle1)
        maty = roty(angle2)
        matz = rotz(angle3)
        #reverse_base_mat = rotz(-base_angle)

        #matrix = np.matmul(maty, matz)
        #matrix = np.matmul(matx, matrix)

        fitting_copy = fitting_model.copy()

        for i in range(fitting_copy.n_out_ports):
            point_delta = fitting_copy.ports_shift[i]
            point_vector = np.array(point_delta)
            #rotated_vector = np.matmul(matrix, point_vector)
            rotated_vector = np.matmul(matz, point_vector)
            rotated_vector = np.matmul(maty, rotated_vector)
            rotated_vector = np.matmul(matx, rotated_vector)
            if direction == 'left':
                if wall.plane.is_perpendicular(Line((0, 0, 0), (1, 0, 0))):
                    rotated_vector[1] *= -1
                elif wall.plane.is_perpendicular(Line((0, 0, 0), (0, 1, 0))):
                    rotated_vector[0] *= -1

            rotated_vector = np.matmul(base_mat, rotated_vector)

            fitting_copy.ports_shift[i] = Point3D(rotated_vector)

        fitting_copy.in_port = in_port

        line = Line(wall.l_point, wall.l_point + Point3D(0, 0, wall.height))

        for i in range(fitting_model.n_out_ports):
            port_copy = in_port.copy()
            port_copy.owner = fitting_copy
            port_copy.diam = fitting_copy.diams[1 + i]
            port_copy.point += fitting_copy.ports_shift[i]

            projection = wall.plane.projection(port_copy.point)
            flat_wall_x = line.distance(projection)
            port_copy.flat_point = Point(flat_wall_x + wall.delta_l, projection.z)

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

    def export_data(self, filename):
        """
        Exports information about pipes, fittings, it's connection
        :return:
        """
        type = []
        id_table = []
        fit_type = []
        fit_id_table = []
        diams = []
        x = []  # in
        y = []  # in
        z = []  # in
        x_fit = []
        y_fit = []
        z_fit = []
        x_out = [[], [], []]
        y_out = [[], [], []]
        z_out = [[], [], []]
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
            type.append(f'Фиттинг {self.fittings[j].type}  {self.fittings[j].id}')
            id_table.append(N + j + 1)
            diams.append(' '.join([str(elem) for elem in self.fittings[j].diams]))
            x.append(round(self.fittings[j].in_port.point.x, 3))
            y.append(round(self.fittings[j].in_port.point.y, 3))
            z.append(round(self.fittings[j].in_port.point.z, 3))

            x_fit.append(round(self.fittings[j].in_port.point.x, 3))
            y_fit.append(round(self.fittings[j].in_port.point.y, 3))
            z_fit.append(round(self.fittings[j].in_port.point.z, 3))
            fit_type.append(f'Фиттинг {self.fittings[j].type}  {self.fittings[j].id}')
            fit_id_table.append(N + j + 1)

            for n_out in range(0,3):
                if len(self.fittings[j].out_ports) > n_out:
                    x_out[n_out].append(round(self.fittings[j].out_ports[n_out].point.x, 2))
                    y_out[n_out].append(round(self.fittings[j].out_ports[n_out].point.y,2))
                    z_out[n_out].append(round(self.fittings[j].out_ports[n_out].point.z,2))
                else:
                    x_out[n_out].append(None)
                    y_out[n_out].append(None)
                    z_out[n_out].append(None)

            self.fittings[j].id_table = id_table[-1]
            if self.fittings[j].reduction is not None:
                type.append(self.fittings[j].reduction.type)
                N += 1
                id_table.append(N + j + 1)
                diams.append(self.fittings[j].reduction.diams)
                x.append(x[-1])
                y.append(y[-1])
                z.append(z[-1])
                x1.append(None)
                x2.append(None)
                y1.append(None)
                y2.append(None)
                z1.append(None)
                z2.append(None)

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

        data2 = {'Звено': connections, 'Диаметр': connections_diam, 'Длина [м]': lengths, 'X1': x1, 'Y1': y1, 'Z1': z1, 'X2': x2,'Y2': y2, 'Z2': z2}

        data3 = {'ИД': fit_id_table, 'Тип': fit_type, 'Вход. X': x_fit, 'Вход. Y': y_fit, 'Вход. Z': z_fit,
                 'Выход0. X': x_out[0], 'Выход0. Y': y_out[0], 'Выход0. Z': z_out[0],
                 'Выход1. X': x_out[1], 'Выход1. Y': y_out[1], 'Выход1. Z': z_out[1],
                 'Выход2. X': x_out[2], 'Выход2. Y': y_out[2], 'Выход2. Z': z_out[2]}
        df1 = pd.DataFrame(data1)
        df2 = pd.DataFrame(data2)
        df3 = pd.DataFrame(data3)

        with pd.ExcelWriter(filename) as writer:
            df1.to_excel(writer, sheet_name='Узлы')
            df2.to_excel(writer, sheet_name='Трубы')
            df3.to_excel(writer, sheet_name='Фиттинги')
        return data1, data2, data3

    def calculate_tubes(self):
        """
        Main function to run calculations
        :return:
        """

        units = self.plumbing_units + [self.water_raiser]

        for unit in units:
            unit.room = self

        corners = self.process_walls(units)

        l_port, l_branch, l_corners, r_port, r_branch, r_corners = self.process_water_raiser(corners)

        self.process_plumbing_branch(r_port, r_branch, r_corners, 'right')
        self.process_plumbing_branch(l_port, l_branch, l_corners, 'left')

        #self.export_data()


#
#
# fittings_allow = {'Triple': [Fitting('Triple', 101, [110, 110, 110], 87,
#                                      [Point3D(0, 0, 0.255), Point3D(0.121, 0, 0.116)]),
#                              Fitting('Triple', 102, [110, 50, 110], 87,
#                                      [Point3D(0, 0, 0.176), Point3D(0.107, 0, 0.088)]),
#                              Fitting('Triple', 103, [50, 50, 50], 87,
#                                      [Point3D(0, 0, 0.15), Point3D(0.077, 0, 0.076)]),
#                              Fitting('Triple', 104, [110, 110, 110], 45,
#                                      [Point3D(0, 0, 0.280), Point3D(0.137, 0, 0.228)]),
#                              Fitting('Triple', 106, [50, 50, 50], 45,
#                                      [Point3D(0, 0, 0.17), Point3D(0.076, 0, 0.14)])
#                              ],
#
#                   'Quadruple': [Fitting('Quadruple', 201, [110, 110, 110, 50], 87,
#                                         [Point3D(0, 0, 0.244), Point3D(0.107, 0, 0.122), Point3D(-0.127, 0, 0.124)])],
#
#                   'Turning': [Fitting('Turning', 310, [40, 40], 87, [Point3D(0.075, 0, 0.08)]),
#                               Fitting('Turning', 305, [50, 50], 87, [Point3D(0.075, 0, 0.08)]),
#                               Fitting('Turning', 306, [110, 110], 87, [Point3D(0.117, 0, 0.123)]),
#                               Fitting('Turning', 309, [40, 40], 45, [Point3D(0.044, 0, 0.107)]),
#                               Fitting('Turning', 304, [50, 50], 45, [Point3D(0.044, 0, 0.107)]),
#                               Fitting('Turning', 300, [110, 110], 45, [Point3D(0.062, 0, 0.150)])],
#
#                   'Reduction': [Fitting('Reduction', 401, [110, 50], 0, [Point3D(0.027, 0, 0.079)]),
#                                 Fitting('Reduction', 402, [110, 40], 0, [Point3D(0.027, 0, 0.079)]),
#                                 Fitting('Reduction', 403, [50, 40], 0, [Point3D(0.002, 0, 0.015)])]}
#
# tubes_allow = {'Tube': [Tube('Tube', 500, 110, 3, 0.06),
#                         Tube('Tube', 501, 50, 3, 0.044),
#                         Tube('Tube', 502, 40, 2, 0.039)],
#                'Nozzle': []}
#
#
# tubes_allow = {}
# fittings_allow = {}


def import_tables(database, room):
    fittings_allow = {'Triple': [Fitting('Triple', 101, [110, 110, 110], 87,
                                         [Point3D(0, 0, 0.255), Point3D(0.121, 0, 0.116)]),
                                 Fitting('Triple', 102, [110, 50, 110], 87,
                                         [Point3D(0, 0, 0.176), Point3D(0.107, 0, 0.088)]),
                                 Fitting('Triple', 103, [50, 50, 50], 87,
                                         [Point3D(0, 0, 0.15), Point3D(0.077, 0, 0.076)]),
                                 Fitting('Triple', 104, [110, 110, 110], 45,
                                         [Point3D(0, 0, 0.280), Point3D(0.137, 0, 0.228)]),
                                 Fitting('Triple', 106, [50, 50, 50], 45,
                                         [Point3D(0, 0, 0.17), Point3D(0.076, 0, 0.14)])
                                 ],

                      'Quadruple': [Fitting('Quadruple', 201, [110, 110, 110, 50], 87,
                                            [Point3D(0, 0, 0.244), Point3D(0.107, 0, 0.122), Point3D(-0.127, 0, 0.124)])],

                      'Turning': [Fitting('Turning', 310, [40, 40], 87, [Point3D(0.075, 0, 0.08)]),
                                  Fitting('Turning', 305, [50, 50], 87, [Point3D(0.075, 0, 0.08)]),
                                  Fitting('Turning', 306, [110, 110], 87, [Point3D(0.117, 0, 0.123)]),
                                  Fitting('Turning', 309, [40, 40], 45, [Point3D(0.044, 0, 0.107)]),
                                  Fitting('Turning', 304, [50, 50], 45, [Point3D(0.044, 0, 0.107)]),
                                  Fitting('Turning', 300, [110, 110], 45, [Point3D(0.062, 0, 0.150)])],

                      'Reduction': [Fitting('Reduction', 401, [110, 50], 0, [Point3D(0.027, 0, 0.079)]),
                                    Fitting('Reduction', 402, [110, 40], 0, [Point3D(0.027, 0, 0.079)]),
                                    Fitting('Reduction', 403, [50, 40], 0, [Point3D(0.002, 0, 0.015)])]}

    tubes_allow = {'Tube': [Tube('Tube', 500, 110, 3, 0.06),
                            Tube('Tube', 501, 50, 3, 0.044),
                            Tube('Tube', 502, 40, 2, 0.039)],
                   'Nozzle': []}
    room.tubes_allow = tubes_allow
    room.fittings_allow = fittings_allow
    # con = sqlite3.connect(database)
    # cur = con.cursor()
    #
    # cur.execute("SELECT * FROM main_material")
    #
    # rows = cur.fetchall()
    #
    # for row in rows:
    #     type = row[0].capitalize()
    #     id = float(row[1])
    #     if type == 'Quadruple':
    #         N = 5
    #     elif type == 'Triple':
    #         N = 4
    #     elif type == 'Turning' or type == 'Reduction':
    #         N = 3
    #     else:
    #         print('Unknown fitting type ', type)
    #         break
    #
    #     diams = [float(row[2])]
    #     for n in range(2, N):
    #         diams.append(float(row[n]))
    #     angle = float(row[5])
    #     ports_shift = []
    #     for i in range(len(diams)-1):
    #         ports_shift.append(Point3D(float(row[6 + 3 * i]), float(row[6 + 3 * i + 1]), float(row[6 + 3 * i + 2])))
    #     room.fittings_allow[type].append(Fitting(type, id, diams, angle, ports_shift))
    #
    # cur.execute("SELECT * FROM tube")
    #
    # rows = cur.fetchall()
    #
    # for row in rows:
    #     id = row[0]
    #     diam = float(row[1])
    #     nominal_length = float(row[2])
    #     connection_shift = float(row[3])
    #     room.tubes_allow['Tube'].append(Tube('Tube', id, diam, nominal_length, connection_shift))
    # con.close()



def initialize_room(plumbing_units, walls):
    room = Room()

    for pu in plumbing_units:
        room.plumbing_units.append(PlumbingUnit(*pu))

    room.water_raiser = PlumbingUnit(0,0,0,110)

    for wall in walls:
        room.add_wall(*wall, 1)
    # room1

    # room.add_wall(-0.15, 1.568, -0.1535, -0.1535, 0.353, 1.2)
    # room.add_wall(1.568, 1.568, -0.1535, -2.310, 0.1, 2)
    # room.add_wall(1.568, -0.15, -2.310, -2.310, 0.05, 2)
    # room.add_wall(-0.15, -0.15, -0.1535, -1, 0.1, 2)
    # #
    # raiser = PlumbingUnit(0, 0, 0, 110)
    # toilet = PlumbingUnit(0.34745, -0.153, 0.2485, 110)
    # sink = PlumbingUnit(1.567, -1.2301, 0.616, 50)
    # bath = PlumbingUnit(1.567, -1.943, 1.19850, 50)
    # boss = PlumbingUnit(1.585069, -2.648512, 0.410237, 40)
    # sink2 = PlumbingUnit(-0.14, -0.5, 0.5, 40)
    # #
    # room.water_raiser = raiser
    # room.plumbing_units = [toilet, sink, bath]
    
    room.tubes_allow = {'Tube': [], 'Nozzle': []}
    room.fittings_allow = {'Triple': [], 'Quadruple': [], 'Turning': [], 'Reduction': []}

    # # room2 cursed
    # raiser = PlumbingUnit(0, 0, 0, 110)
    # toilet = PlumbingUnit(0.591, -0.154, 0.236, 110)
    # sink = PlumbingUnit(1.243, -1.236, 0.753, 50)
    #
    # bath1 = PlumbingUnit(-1.9, -0.6, 0.4, 50)
    # bath2 = PlumbingUnit(-1.9, -1, 0.2, 50)
    #
    # hole = PlumbingUnit(-1.2, -1.5, 0.6, 50)
    # hole.direct = True
    #
    # room.add_wall(-1.193, 1.221, -0.1765, -0.1765, 0.292, 0.290)
    # room.add_wall(1.221,  1.221, -0.1765, -1.686, 0.072, 2)
    # room.add_wall(-1.193,  -1.193, -1.686, -0.1765,  0.072, 2)
    #
    # room.water_raiser = raiser
    # room.plumbing_units = [toilet, sink, bath1, bath2, hole]

    # # room3
    # room.add_wall(-0.15, 1.568, -0.1535, -0.1535, 0.353, 1.2)
    # room.add_wall(1.568, 1.568, -0.1535, -2.310, 0.1, 2)
    # room.add_wall(-0.15, -0.15, -0.1535, -2.31, 0.1, 2)
    #
    # room.add_wall(1.568, 3.52, -0.1535, -0.1535, 0.1, 2)
    # room.add_wall(3.52, 3.52, -0.1535, -2.31, 0.1, 2)
    #
    # raiser = PlumbingUnit(0, 0, 0, 110)
    # toilet = PlumbingUnit(0.34745, -0.153, 0.2485, 110)
    #
    # sink = PlumbingUnit(1, -0.153, 0.616, 50)
    # bath = PlumbingUnit(3.51, -1.943, 1.19850, 50)
    #
    # #boss = PlumbingUnit(1.585069, -2.6409, 0.410237, 40)
    # sink2 = PlumbingUnit(1.57, -1.5, 0.5, 50)
    # sink3 = PlumbingUnit(-3.51, -1.5, 0.8, 50)
    #
    # room.double = True
    # #boss.direct = True
    #
    # room.water_raiser = raiser
    # room.plumbing_units = [toilet, sink, bath, sink2, sink3]

    room.tubes_allow = {'Tube': [], 'Nozzle': []}
    room.fittings_allow = {'Triple': [], 'Quadruple': [], 'Turning': [], 'Reduction': []}

    return room

#


walls = [(-0.15, 1.568, -0.1535, -0.1535, 0.353), (1.568, 1.568, -0.1535, -2.310, 0.1),
         (1.568, -0.15, -2.310, -2.310, 0.05), (-0.15, -0.15, -0.1535, -1, 0.1)]
PU = [(0.34745, -0.153, 0.2485, 110), (1.567, -1.2301, 0.616, 50), (1.567, -1.943, 0.35850, 50)]

room1 = initialize_room(PU, walls)

import_tables('db.sqlite3', room1)

room1.calculate_tubes()

room1.export_data('dataframe2.xlsx')

