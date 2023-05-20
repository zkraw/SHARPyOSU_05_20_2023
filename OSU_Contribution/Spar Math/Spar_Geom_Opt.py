import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from shapely.geometry import Polygon
from scipy.optimize import minimize
import sys

# NACA63 = reader('NACA63015F.dat')
class Spar_Max():
    def __init__(self, Name, Front = 0, Back = 1):
        self.airfoil = self.reader(Name)
        self.poly = Polygon(self.airfoil)
        self.start = Front
        self.end = Back

        if self.start < 0 or self.start > 1:
            print('ERROR: Enter front bound as a percentage')
            sys.exit()
        elif self.end < 0 or self.end > 1:
            print('ERROR: Enter back bound as a percentage')
            sys.exit()
        elif self.start == 1 or self.end == 0 or self.start >= self.end:
            print('ERROR:Enter in valid front and back bound combo')
            sys.exit()

    def reader(self, Name):
        airfoil = []
        f = open(Name, "r")
        for line in f:
            current = line.strip('\n')
            current = current.split(" ")
            final = []
            for i in range(len(current)):
                if current[i] != '':
                    final.append(current[i])
            current = final
            try:
                one = float(current[0])
                two = float(current[1])
                if type(float(current[1])) == float:
                    airfoil.append([one, two])
            except:
                pass
        return np.array(airfoil)

    def interp(self, x, y, middle):
        answer = (y[1] - y[0]) / (x[1] - x[0]) * (middle - x[0]) + y[0]
        return answer

    def y_boundary(self, x):
        x_ref = self.poly.exterior.xy[0]
        y = self.poly.exterior.xy[1]
        want = []

        check = [0, 0]
        for i in range(len(x_ref)):
            if x == x_ref[i]:
                want.append(y[i])
                if len(want) > 2:
                    want[1] = y[i]
                    break
            elif x > x_ref[i]:
                min = x_ref[i]
                low = y[i]
                check[0] = 1
                if check == [1, 1]:
                    want.append(self.interp([min, max], [low, high], x))
                    check = [0, 0]

            elif x < x_ref[i]:
                max = x_ref[i]
                high = y[i]
                check[1] = 1
                if check == [1, 1]:
                    want.append(self.interp([min, max], [low, high], x))
                    check = [0, 0]
                 # thick in mm

        return want

    def minimize_circle(self):
        poly = self.poly
        half_x = (poly.bounds[2] + poly.bounds[1]) / 2
        self.theta_high = np.linspace(0, np.pi, 15)
        self.theta_low = np.linspace(np.pi, 2*np.pi, 15)
        self.x0 = half_x
        self.r = .01

        def area(x):
            return -np.pi*x[2]**2

        class constraint(Spar_Max):
            def __init__(self, theta, poly):
                self.theta = theta
                self.poly = poly
            def polar_to_cart(self, origin, theta, r):
                x = origin[0] + r*np.cos(theta)
                y = origin[1] + r*np.sin(theta)
                return x, y

            def constraint_below(self, x):
                xx, yy = self.polar_to_cart((x[0], x[1]), self.theta, x[2])
                if xx <= 0:
                    xx = 0
                elif xx >= 1:
                    xx = 1
                y_bound = min(self.y_boundary(xx))
                return yy - y_bound

            def constraint_above(self, x):
                xx, yy = self.polar_to_cart((x[0], x[1]), self.theta, x[2])
                if xx <= 0:
                    xx = 0
                elif xx >= 1:
                    xx = 1
                y_bound = max(self.y_boundary(xx))
                return y_bound - yy


        cons = []
        for jj in range(len(self.theta_low)):
            cons.append({'type': 'ineq', 'fun': constraint(self.theta_low[jj], self.poly).constraint_below})
        for jj in range(len(self.theta_high)):
            cons.append({'type': 'ineq', 'fun': constraint(self.theta_high[jj], self.poly).constraint_above})

        if self.start == 0 and self.end == 1:
            self.bnds = ((poly.bounds[0], poly.bounds[2]),
                    (poly.bounds[1], poly.bounds[3]),
                    (0, poly.bounds[2]))

        elif self.start != 0 and self.end == 1:
            distance = (poly.bounds[2] - poly.bounds[0]) * self.start
            self.bnds = ((distance + self.r, poly.bounds[2] -self.r),
                    (poly.bounds[1], poly.bounds[3]),
                    (0, min((self.x0 - poly.bounds[0], poly.bounds[2] - self.x0))))

        elif self.start == 0 and self.end != 1:
            end = (poly.bounds[2] - poly.bounds[0]) * self.end
            self.bnds = ((poly.bounds[0] + self.r, end - self.r),
                    (poly.bounds[1], poly.bounds[3]),
                    (0, min((self.x0 - poly.bounds[0], poly.bounds[2] - self.x0))))

        elif self.start != 0 and self.end != 1:
            front = (poly.bounds[2] - poly.bounds[0]) * self.start
            end = (poly.bounds[2] - poly.bounds[0]) * self.end
            self.bnds = ((front + self.r, end - self.r),
                    (poly.bounds[1], poly.bounds[3]),
                    (0, min((self.x0 - poly.bounds[0], poly.bounds[2] - self.x0))))

        # [center x, center y, radius]
        half_y = (poly.bounds[3] + poly.bounds[1]) / 2
        x0 = [half_x, half_y, .2]
        self.results = minimize(area, x0, constraints=cons,
                           method='SLSQP', bounds= self.bnds)

        theta = np.linspace(0, 2*np.pi, 100)
        final = np.zeros((len(theta), 2))
        x0 = self.results.x[0]
        y0 = self.results.x[1]
        r = self.results.x[2]

        for i in range(len(theta)):
            final[i, :] = [r*np.cos(theta[i]) + x0, r*np.sin(theta[i]) + y0]

        circle = Polygon(final)

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 2))
        ax.plot(*poly.exterior.xy, 'r')
        ax.plot(*circle.exterior.xy, 'b')
        plt.fill_between(*circle.exterior.xy)
        plt.show()

    def minimize_rectangle(self):
        # points = (x0, x1, y0, y1)
        def area(x):
            x1 = x[1]
            x0 = x[0]
            y1 = x[3]
            y0 = x[2]
            area = abs((x1 - x0) * (y1 - y0))
            return -area

        class constraint(Spar_Max):
            def __init__(self, percent, poly):
                self.percent = percent
                self.poly = poly

            def bottom(self, x):
                per = self.percent
                xx = x[1] - x[0]
                yy = x[2]
                x_temp = per*xx + x[0]
                y_bound = min(self.y_boundary(x_temp))
                return yy - y_bound



        # constraints are under here
        # NEED to add another constraint that checks if all points on horizontal line is
        # above and below bounds
        # constraints for first corner much be within the airfoil (lb < y0 <= ub )
        def constraint1(x):
            hold = self.y_boundary(x[0])
            random = min(hold)
            return x[2] - random

        def constraint2(x):
            hold = self.y_boundary(x[0])
            random = max(hold)
            return random - x[2]

        # constraints for second corner much be within the airfoil (lb < y0 <= ub )
        def constraint3(x):
            hold = self.y_boundary(x[1])
            random = min(hold)
            return x[3] - random

        def constraint4(x):
            hold = self.y_boundary(x[1])
            random = max(hold)
            return random - x[3]

        # third corner which is above first
        def constraint5(x):
            hold = self.y_boundary(x[0])
            random = min(hold)
            return x[3] - random

        def constraint6(x):
            hold = self.y_boundary(x[0])
            random = max(hold)
            return random - x[3]

        # fourth corner which is below second
        def constraint7(x):
            hold = self.y_boundary(x[1])
            random = min(hold)
            return x[2] - random

        def constraint8(x):
            hold = self.y_boundary(x[1])
            random = max(hold)
            return random - x[2]

        con1 = {'type': 'ineq', 'fun': constraint1}
        con2 = {'type': 'ineq', 'fun': constraint2}
        con3 = {'type': 'ineq', 'fun': constraint3}
        con4 = {'type': 'ineq', 'fun': constraint4}
        con5 = {'type': 'ineq', 'fun': constraint5}
        con6 = {'type': 'ineq', 'fun': constraint6}
        con7 = {'type': 'ineq', 'fun': constraint7}
        con8 = {'type': 'ineq', 'fun': constraint8}
        cons = [con1, con2, con3, con4, con5, con6, con7, con8]
        per = np.linspace(.01, .99, 10)
        for i in range(len(per)):
            cons.append({'type': 'ineq', 'fun': constraint(per[i], self.poly).bottom})


        poly = self.poly
        half = (poly.bounds[3] + poly.bounds[1]) / 2
        if self.start == 0 and self.end == 1:
            bnds = ((poly.bounds[0], poly.bounds[2]), (poly.bounds[0], poly.bounds[2]),
                         (poly.bounds[1], half), (half, poly.bounds[3]))
        elif self.start != 0 and self.end == 1:
            distance = (poly.bounds[2] - poly.bounds[0])* self.start
            bnds = ((distance, poly.bounds[2]), (poly.bounds[0], poly.bounds[2]),
                    (poly.bounds[1], half), (half, poly.bounds[3]))

        elif self.start == 0 and self.end != 1:
            end = (poly.bounds[2] - poly.bounds[0])* self.end
            bnds = ((poly.bounds[0], poly.bounds[2]), (poly.bounds[0], end),
                    (poly.bounds[1], half), (half, poly.bounds[3]))

        elif self.start != 0 and self.end != 1:
            front = (poly.bounds[2] - poly.bounds[0])* self.start
            end = (poly.bounds[2] - poly.bounds[0])* self.end
            bnds = ((front, poly.bounds[2]), (poly.bounds[0], end), (poly.bounds[1], half), (half, poly.bounds[3]))

        x0 = [0.1, .99, 0, 1]
        results = minimize(area, x0, constraints= cons,
                                      method = 'SLSQP', bounds = bnds)

        final = np.zeros((4,2))
        final[0, :] = [results.x[0], results.x[2]]
        final[1, :] = [results.x[0], results.x[2] + (results.x[3] - results.x[2])]
        final[2, :] = [results.x[1], results.x[3]]
        final[3, :] = [results.x[1], results.x[3] - (results.x[3] - results.x[2])]

        rect = Polygon(final)

        # fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 2))
        # ax.plot(*poly.exterior.xy, 'r')
        # ax.plot(*rect.exterior.xy, 'b')
        # plt.fill_between(*rect.exterior.xy)
        # plt.show()
        return poly, rect
    
    def minimize_rectangle_center(self, center):
        if center > 1 or center < 0:
            print('ERROR: Center location must be between 0 and 1')
            sys.exit()
        
        global center_loc
        center_loc = center
        def area(x):
            xc = center_loc
            yc = x[0]
            x1 = x[1]
            y1 = x[2]
            area = 2*(y1 - yc) * 2*(x1 - xc)
            return -area

            
        # x1 must be greater than the center location
        def constraint1(x):
            xc = center_loc
            x1 = x[1]
            return x1 - xc
        
        # y1 must be greater than the center location
        def constraint2(x):
            yc = x[0] 
            y1 = x[2]
            return y1 - yc
        
        # yc must be in airfoil 
        def constraint3(x):
            yc = x[0]
            hold = self.y_boundary(center_loc)
            random = min(hold)
            return yc - random

        def constraint4(x):
            yc = x[0]
            hold = self.y_boundary(center_loc)
            random = max(hold)
            return random - yc 
            
        # y1 must be below upper airfoil bound 
        def constraint5(x):
            x1 = x[1]
            y1 = x[2]
            hold = self.y_boundary(x1)
            random = max(hold)
            return random - y1 
       
        # point below corner must be above airfoil lower bound 
        def constraint6(x):
            yc = x[0]
            x1 = x[1]
            y1 = x[2]
            hold = self.y_boundary(x1)
            random = min(hold)
            height = 2 * (y1 - yc)
            return (y1 - height) - random

        # upper left point must be below airfoil upper bound 
        def constraint7(x):
            x1 = x[1]
            y1 = x[2]
            base = 2*(x1 - center_loc)
            x2 = x1 - base
            hold = self.y_boundary(x2)
            random = max(hold)
            return random - y1 
        
        # lower left point must be above airfoil lower bound 
        def constraint8(x):
            yc = x[0]
            x1 = x[1]
            y1 = x[2]
            base = 2*(x1 - center_loc)
            height = 2*(y1 - yc)
            x2 = x1 - base
            y2 = y1 - height
            hold = self.y_boundary(x2)
            random = min(hold)
            return y2 - random 

        class constraint(Spar_Max):
            def __init__(self, percent, poly):
                self.percent = percent
                self.poly = poly

            def bottom(self, x):
                per = self.percent
                xc = center_loc
                yc = x[0]
                x1 = x[1] 
                y1 = x[2]
                base = 2 * (x1 - xc)
                height = 2 * (y1 - yc)
                x_temp = x1 - per*base
                y_bound = min(self.y_boundary(x_temp))
                return (yc - height) - y_bound

            def top(self, x):
                per = self.percent
                xc = center_loc
                x1 = x[1] 
                y1 = x[2]
                base = 2 * (x1 - xc)
                x_temp = x1 - per*base
                y_bound = min(self.y_boundary(x_temp))
                return y_bound - y1

        con1 = {'type': 'ineq', 'fun': constraint1}
        con2 = {'type': 'ineq', 'fun': constraint2}
        con3 = {'type': 'ineq', 'fun': constraint3}
        con4 = {'type': 'ineq', 'fun': constraint4}
        con5 = {'type': 'ineq', 'fun': constraint5}
        con6 = {'type': 'ineq', 'fun': constraint6}
        con7 = {'type': 'ineq', 'fun': constraint7}
        con8 = {'type': 'ineq', 'fun': constraint8}
        cons = [con1, con2, con3, con4, con5, con6, con7, con8]
        
        per = np.linspace(.01, .99, 10)
        # for i in range(len(per)):
        #     cons.append({'type': 'ineq', 'fun': constraint(per[i], self.poly).bottom})
        #     cons.append({'type': 'ineq', 'fun': constraint(per[i], self.poly).top})

        diff = 1 - center_loc
        if center_loc <= diff:
            x_bound = center_loc * 2
        else:
            x_bound = self.poly.bounds[2]

        bnds = ((self.poly.bounds[1], self.poly.bounds[3]), (center_loc, self.end), (self.poly.bounds[1], self.poly.bounds[3]))
        x0 = [0, .5, .1]
        results = minimize(area, x0, constraints= cons,
                                      method = 'SLSQP', bounds = bnds)

        final = np.zeros((4,2))
        tolerance = .9
        base = 2 * (results.x[1] - center_loc)
        height = 2 * (results.x[2] - results.x[0])
        
        # final four points
        final[0, :] = [results.x[1], results.x[2]- (1-tolerance)*.5 * height]
        final[1, :] = [results.x[1], results.x[2] - height + (1-tolerance)*.5 * height]
        final[2, :] = [results.x[1] - base, results.x[2] - height + (1-tolerance)*.5 * height]
        final[3, :] = [results.x[1] - base, results.x[2] - (1-tolerance)*.5 * height]

        miny = np.zeros((4,2))
        miny[0, :] = [final[0][0] - .0025, final[0][1] - .0025]
        miny[1, :] = [final[1][0] - .0025, final[1][1] + .0025]
        miny[2, :] = [final[2][0] + .0025, final[2][1] + .0025]
        miny[3, :] = [final[3][0] + .0025, final[3][1] - .0025]

        rect = Polygon(final)
        inner = Polygon(miny)

        return self.poly, rect, inner
         


# poly, rect = Spar_Max('/mnt/c/Users/Zack/Desktop/python/sharpy/Zack_Beam/Spar Math/naca0015.dat', .2, .36).minimize_rectangle_center(.3)
poly, rect, inner = Spar_Max('C:/Users/Zack/Desktop/python/sharpy/Zack_Beam/Spar Math/naca0015.dat', .2, .36).minimize_rectangle_center(.3)


print('Base / chord = ' + str(rect.bounds[2] - rect.bounds[0]))
print('\nHeight / chord = ' + str(rect.bounds[3] - rect.bounds[1]))

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 2))
ax.plot(*poly.exterior.xy, 'r')
ax.plot(*rect.exterior.xy, 'black')
ax.plot(*inner.exterior.xy, 'black')
ax.scatter(.3, 0, 50, color = 'black', marker = 'x')
ax.tick_params(axis= 'x', labelsize = 12)
ax.tick_params(axis= 'y', labelsize = 12)
# plt.fill_between(*rect.exterior.xy)
plt.xlabel('Unit Chord Length', fontsize = 14)
plt.ylabel('Unit Height', fontsize = 14)
fig.tight_layout()
plt.show()


