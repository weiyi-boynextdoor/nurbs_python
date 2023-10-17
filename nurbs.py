import matplotlib.pyplot as plt
import copy

class Nurbs:
    def __init__(self, control_point_num, order) -> None:
        self.control_point_num = control_point_num
        self.order = order
        self.degree = order - 1
        self.cyclic = False
        self.weights = [1 for i in range(control_point_num)]
        # default knot vector
        self.knot_vector = [i for i in range(control_point_num + order)]

    def set_weights(self, weights : list):
        assert(len(weights) == self.control_point_num)
        self.weights = copy.deepcopy(weights)

    def get_basis_value(self, t : float):
        start = 0
        end = 0
        for i in range(self.control_point_num + self.degree):
            knots_equal = self.knot_vector[i] == self.knot_vector[i+1]
            if knots_equal or t < self.knot_vector[i] or t > self.knot_vector[i+1]:
                continue
            start = max(i - self.degree, 0)
            end = i
            break

        buffer = [0.0 for i in range(self.order * 2)]
        buffer[end - start] = 1.0

        for i_order in range(2, self.order + 1):
            if end + i_order >= len(self.knot_vector):
                end = self.control_point_num + self.degree - i_order
            for i in range(end -start + 1):
                knot_index = start + i
                new_basis = 0.0

                if self.knot_vector[knot_index + i_order - 1] - self.knot_vector[knot_index] == 0:
                    assert(buffer[i] == 0)
                if  self.knot_vector[knot_index + i_order] - self.knot_vector[knot_index + 1] == 0:
                    assert(buffer[i+1] == 0)
                
                if buffer[i] != 0.0:
                    new_basis += (t - self.knot_vector[knot_index]) * buffer[i] / (self.knot_vector[knot_index + i_order - 1] - self.knot_vector[knot_index])
                if buffer[i+1] != 0.0:
                    new_basis += (self.knot_vector[knot_index + i_order] - t) * buffer[i+1] / (self.knot_vector[knot_index + i_order] - self.knot_vector[knot_index + 1])
                buffer[i] = new_basis

        return {"start": start, "values": buffer[0:self.order]}

    def get_overall_t_range(self):
        return [self.knot_vector[0], self.knot_vector[-1]]

    def get_valid_t_range(self):
        return [self.knot_vector[self.degree], self.knot_vector[-self.order]]

    def set_template(self, is_end_point, is_bezier):
        repeat_inner = self.degree if is_bezier else 1
        i = 0
        t = 0.0
        if is_end_point:
            # use self.order 0's as header
            while i < self.order:
                self.knot_vector[i] = t
                i += 1
        else:
            self.knot_vector[0] = t
            i += 1
        
        tail_len = self.order if is_end_point else 0
        while i < len(self.knot_vector) - tail_len:
            t += 1.0
            for j in range(repeat_inner):
                self.knot_vector[i] = t
                i += 1

        t += 1.0
        while i < len(self.knot_vector):
            self.knot_vector[i] = t
            i += 1

    def set_knot_vector(self, knot_vector):
        assert(len(knot_vector) == len(self.knot_vector))
        self.knot_vector = copy.deepcopy(knot_vector)

if __name__ == "__main__":
    n_points = 5
    order = 3

    nurbs = Nurbs(n_points, order)
    nurbs.set_template(False, False)
    # nurbs.set_knot_vector([0,1,2,3.5,4,5,6,7,8])
    
    overall_t_range = nurbs.get_overall_t_range()
    valid_t_range = nurbs.get_valid_t_range()

    t_range = overall_t_range

    samples = 401
    step = float(t_range[1] - t_range[0]) / (samples - 1)

    x = [t_range[0] + step * i for i in range(samples)]
    x[-1] = t_range[1]
    
    y_s = []
    for i in range(n_points):
        y_s.append([0.0 for j in range(samples)])

    for i_x in range(len(x)):
        t = x[i_x]
        start_and_values = nurbs.get_basis_value(t)
        start = start_and_values["start"]
        values = start_and_values["values"]
        for i_v in range(len(values)):
            if start + i_v >= n_points:
                continue
            y_s[start + i_v][i_x] = values[i_v]

    plt.figure(figsize=(12,6))
    for i in range(n_points):
        plt.plot(x, y_s[i], label="B%d"%i)
    plt.legend()    
    plt.show()