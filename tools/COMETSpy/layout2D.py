from __future__ import division
import math
sys.path.append('../../')
from tools.customError import *
from meshPoint import *

class layout2D(object):
    """
    Defines a particular two-dimensional layout for the simulation of heterogenous condition 

    Ali R. Zomorrodi - Boston University
    Last updated: July 02 2014
    """

    def __init__(self,layoutType, Lx = None, Ly = None, R = None, h = None):
       
        """
        INPUTS:
        ------
          layoutType: Type of the layout. Eligible cases are rectangle and circle
        boundaryCond: Boundary condition. A dictionary whose keys show the 
                      boundary (different on the layout) and values are tuples 
                      where the first element is a string representing the type 
                      of boundary condition and the second element shows its value.
                      Types of boundary condition include: Dirichlet and Neumannn.
                      The format of inputing the boundary conditions for each 
                      layout is as follows:

                      Rectangle:
                         {'x=0':['type',value],'x=Ly':['type',value],
                          'y=0':['type',value],'y=Lx':['type',value]}
                         Example: {'x=0':['Drichlet',5],'x=10':['Drichlet',5],
                                   'y=0':['Neumann',0],'y=10':['Neumann',0]}
                      Circle: 
                         {'r=R':['type',value]}  
                         Example: {'r=10':['Neumann',0]}

                  Lx: Length of the rectangle in the x-direction 
                      (Optional input. Default: 1 cm)
                  Lx: Length of the rectangle in the y-direction 
                      (Optional input. Default: 1 cm)
                   R: Redius of the circle (Default: 1 cm)
                   h: Mesh size (Optional input. Default: 0.1 cm)
        """

        # Type of the layout. Eligible cases are rectangle and circle
        if layoutType.lower() in ['rectangle','circle']:  
            self.layoutType = layoutType
        else:
            raise customError("Error in defining the layout! Eligible layouts are 'rectangle' and 'circle'\n\n")

        # Mesh size
        if h == None:
            self.h = 0.1
        else:
            self.h = h

        # Rectangle 
        if self.layoutType.lower() == 'rectangle': 
            # Length of the rectangle in the x-direction
            if Lx == None:
                print 'Lx was set to the default value of 1 cm\n'
                self.Lx = 1
            else:
                self.Lx = Lx

            # Length of the rectangle in the y-direction
            if Ly == None:
                print 'Ly was set to the default value of 1 cm\n'
                self.Ly = 1
            else:
                self.Ly = Ly

        # Circle 
        elif self.layoutType.lower() == 'circle':

            if R == None:
                print 'R was set to the default value of 1 cm\n'
                self.R = 1
            else:
                self.R = R 

    #--- Rectangle ---
    def defineRectangleLayout(self):
        """ 
        This method defines the mesh points for a rectangle, i.e., the position 
        of all mesh points, position of internal mesh ooints and the position of 
        mesh points on the boundary

        The origin of a rectangular cooerdinate is set on the lower left vertex
        """

        #--- Define a 2-D mesh ---
        # Divide x- and y-axis
        self.xPoints = self.frange(0,self.Lx,self.h)
        self.yPoints = self.frange(0,self.Ly,self.h)

        # Position (xy-coordinate) of boundary points
        boundary_xyCoord = [(0,j) for j in self.yPoints] + [(self.Lx,j) for j in self.yPoints] + [(i,0) for i in self.xPoints[1:len(self.xPoints)-1]] + [(i,self.Ly) for i in self.xPoints[1:len(self.xPoints)-1]]  

        # Define the dictionary containing boundary points
        for k in boundary_xyCoord:
            x = k[0]
            y = k[1]
            xLabel = self.xPoints.index(x)
            yLabel = self.yPoints.index(y)
            self.boundaryPoints[(xLabel,yLabel)] = meshPoint(type = 'boundary',x = x, y = y, xLabel = xLabel, yLabel = yLabel) 
                   
        # Position (xy-coordinate) of internal mesh points
        internal_xyCoord = [(i,j) for i in self.xPoints[1:len(self.xPoints)-1] for j in self.yPoints[1:len(self.yPoints)-1]] 

        # Define the dictionary containing internal points
        for k in internal_xyCoord:
            x = k[0]
            y = k[1]
            xLabel = self.xPoints.index(x)
            yLabel = self.yPoints.index(y)
            self.internalPoints[(xLabel,yLabel)] = meshPoint(type = 'internal',x = x, y = y, xLabel = xLabel, yLabel = yLabel) 

        # Now that we have assigned the labels we can define fE, fW, fN and fS
        self.fCalc()


    #--- Circle ---
    def defineCircleLayout(self):
        """ 
        This method defines the mesh points for a circle, i.e., the position 
        of all mesh points, position of internal mesh ooints and the position 
        of mesh points on boundary

        The circle is assumed to be circumvented in a square whose side is 2R,
        where R is the radius of the circle. 
        The origin of coordinate is put on the lower left vertex of this square. 
        Note that if the origin was in the center of circle its equation will 
        be x^2 + y^2 = R^2. If the origin is moved to the lower left vertex its
        equation will be  (x - R)^2 + (y - R)^2 = R^2
        """
        # Define a 2-D array representing the position of each mesh point
        self.xPoints = self.frange(0,self.R,self.h)
        self.yPoints = self.frange(0,self.R,self.h)

        # Position of internal mesh points
        internal_xyCoord = [(i,j) for i in self.xPoints for j in self.yPoints if (i - self.R)**2 + (j - self.R)**2 < self.R^2] 

        # Define the dictionary containing internal points
        for k in internal_xyCoord:
            x = k[0]
            y = k[1]
            xLabel = xPoints.index(x)
            yLabel = yPoints.index(y)
            self.internalPoints[(xLabel,yLabel)] = meshPoint(type = 'internal',x = x, y = y, xLabel = xLabel, yLabel = yLabel) 

        # Position of the boundary points
        # Find the intersection of each mesh line with the circle
        # For a given vertical mesh line: 
        # y = R - sqrt(R^2 - (x-R)^2) & y = R + sqrt(R^2 - (x-R)^2)
        # For a given horizontal mesh line: 
        # x = R - sqrt(R^2 - (y-R)^2) & x = R + sqrt(R^2 - (y-R)^2)
        boundary_xyCoord = [(0,self.R),(self.R,0),(self.R,2*self.R),(2*self.R,self.R)]  + [(x,self.R - math.sqrt(self.R**2 - (x-self.R)**2)) for x in self.xPoints[1:len(self.xPoints)-1]] + [(x,self.R - math.sqrt(self.R**2 + (x-self.R)**2)) for x in self.xPoints[1:len(self.xPoints)-1]] + [(self.R - math.sqrt(self.R**2 - (y-self.R)**2),y) for y in self.yPoints[1:len(yPoints)-1]] + [(self.R + math.sqrt(self.R**2 - (y-self.R)**2),y) for y in self.yPoints[1:len(yPoints)-1]] 

        # Define the dictionary containing boundary points
        for k in boundary_xyCoord:
            x = k[0]
            y = k[1]
            [xLabel,yLabel] = self.findLabel(x,y)
            self.boundaryPoints[(xLabel,yLabel)] = meshPoint(type = 'boundary',x = x, y = y, xLabel = xLabel, yLabel = yLabel) 
                   
        # Now that we have assigned the labels we can define fE, fW, fN and fS
        self.fCalc()

        
    #----- A function for creating a range with fractional increment ------
    def frange(self,start,stop,step):
        """
        Similar to range function of python with two differences. (1) It accepts a 
        fractional increment. (2) the sopt point is always the last element of 
        the list 

        Example frange(0,1.0.3) = [0,0.3,0.6,0.9,1]
        """

        frange_res = []
     
        current_num = start
        while current_num <= stop:
            frange_res.append(current_num)
            current_num = current_num + step

        if frange_res[len(frange_res) - 1] < stop:
            frange_res.append(stop)

        return frange_res

    #----- A function for finding the label of a mesh point ------
    def findLabel(self, x, y):
        """
        Finds the label of a point in the mesh. Here, label simply 
        means numbering of the points.  
        x and y are the x- and y-coordinates of the mesh point 
        """
        #- xLabel -
        if x in self.xPoints:
            xLabel = self.xPoints.index(self.boundary_xyCoord[k][0])
        else:
            # Find the points on the left and right of this point in xPoints
            done = 0
            k = 0
            while done == 0:
                if x > self.xPoints[k]: 
                    done = 1;
                    leftP = (self.xPoints[k],y)
                    rightP = (self.xPoints[k+1],y)
                else:
                    k += 1 

            # Check whether the left or right point is an internal mesh point
            if leftP in self.internal_xyCoord and rightP in self.internal_xyCoord:
                # It is impossible to have both the left and right points inside a 
                # grid unless we deal with a non-convex shape (thid should be 
                # addressed in future developments) 
                raise customError('**Error! Both the left or right points of the given point are inside the mesh') 
            elif leftP in self.internal_xyCoord:
                # If the left point is inside the grid the label of the current i
                # point should be one plus the label of the point on its left
                xLabel = self.xPoints.index(leftP[0]) + 1
            elif rightP in self.internal_xyCoord:
                # If the right point is inside the grid the label of the current 
                # point should be one minos the label of the point on its right
                xLabel = self.xPoints.index(rightP[0]) - 1
            else:
                raise customError('**Error! Neither the left or right points of the given point are inside the mesh') 

        #- yLabel -
        if y in self.self.yPoints:
            yLabel = self.yPoints.index(self.boundary_xyCoord[k][0])
        else:
            # Find the points on the left and right of this point in xPoints
            done = 0
            k = 0
            while done == 0:
                if y > self.yPoints[k]: 
                    done = 1;
                    lowerP = (x,self.yPoints[k])
                    upperP = (x,self.yPoints[k+1])
                else:
                    k += 1 

            # Check whether the lower or upper point is an internal mesh point
            if lowerP in self.internal_xyCoord and upperP in self.internal_xyCoord:
                # It is impossible to have both the lower and upper points 
                # inside a grid unless we deal with a non-convex shape 
                raise customError('**Error! Both the lower or upper points of the given point are inside the mesh') 
            elif lowerP in self.internal_xyCoord:
                # If the lower point is inside the grid the label of the current 
                # point should be one plus the label of the point on its lower
                yLabel = self.xPoints.index(lowerP[0]) + 1
            elif upperP in self.internal_xyCoord:
                # If the upper point is inside the grid the label of the current 
                # point should be one minos the label of the point on its upper
                yLabel = self.xPoints.index(upperP[0]) - 1
            else:
                raise customError('**Error! Neither the lower or upper points of the given point are inside the mesh') 

    #----- A function for finding fE, fW, fN and fS ------
    def fCalc(self):
        """
        Computes fE, fW, fN and fS for each mesh point. These are
        used to calcualte the coefficients needed to computing the Laplacian
        
        INPUTS: self.internalPoints
        OUTPUT: fE, fW, fN and fS are assigned to each instance of the class
                meshPoint stored in self.internalPoints
        """
        # A dictionary composed of all internal and boundary points
        allPoints = dict(self.internalPoints.items() + self.boundaryPoints.items())

        for pointLabel in allPoints.keys():
            # Compute fE, fW, fN and fW only for internal mesh points
            if allPoints[pointLabel].type.lower() == 'internal':
                xLabel = pointLabel[0]
                yLabel = pointLabel[1]
                x = self.internalPoints[(xLabel,yLabel)].x
                y = self.internalPoints[(xLabel,yLabel)].y
                xE = allPoints[(xLabel + 1,yLabel)].x
                xW = allPoints[(xLabel - 1,yLabel)].x
                yN = allPoints[(xLabel,yLabel + 1)].y
                yS = allPoints[(xLabel,yLabel - 1)].y
    
                if (xE - x)/self.h < -0.000001 or (xE - x)/self.h > 1.000001:
                    errorMessage = '**Error! (xE - x)/h for the point with label ('+str(xLabel) + ',' + str(yLabel) + ') and coordinate (' + str(x) + ',' + str(y) + ') is not between zero and one: (xE - x)/h = ' + str((xE - x)/self.h)
                    raise customError(errorMessage)
                else:
                    self.internalPoints[(xLabel,yLabel)].fE = (xE - x)/self.h

                # Note that in the following we use -0.000001 and 1.000001  
                # instead of 0 and 1, respectively, to avoid problems with
                # with very small fractions. For example if the fractions is
                # greater than one by 2.22e-16 the condition (x - xW)/self.h > 1
                # will be false and the code returns an error
                if (x - xW)/self.h < -0.000001 or (x - xW)/self.h > 1.000001:
                    errorMessage = '**Error! (x - xW)/h for the point with label ('+str(xLabel) + ',' + str(yLabel) + ') and coordinate (' + str(x) + ',' + str(y) + ') is not between zero and one: (x - xW)/h = ' + str((x - xW)/self.h)
                    raise customError(errorMessage)
                else:
                    self.internalPoints[(xLabel,yLabel)].fW = (x - xW)/self.h
    
                if (yN - y)/self.h < -0.000001 or (yN - y)/self.h > 1.000001:
                    errorMessage = '**Error! (yN - y)/h for the point with label ('+str(xLabel) + ',' + str(yLabel) + ') and coordinate (' + str(x) + ',' + str(y) + ') is not between zero and one: (yN - y)/h = ' + str((yN - y)/self.h)
                    raise customError(errorMessage)
                else:
                    self.internalPoints[(xLabel,yLabel)].fN = (yN - y)/self.h
    
                if (y - yS)/self.h < -0.000001 or (y - yS)/self.h > 1.000001:
                    errorMessage = '**Error! (y - yS)/h for the point with label ('+str(xLabel) + ',' + str(yLabel) + ') and coordinate (' + str(x) + ',' + str(y) + ') is not between zero and one: (y - yS)/h = ' + str((y - yS)/self.h)
                    raise customError(errorMessage)
                else:
                    self.internalPoints[(xLabel,yLabel)].fS = (y - yS)/self.h
    
                # Calculate the coeeficients requried to compute the Laplacian  
                self.internalPoints[(xLabel,yLabel)].LapCoeffCalc()


    def run(self):
        """ 
        This method defines the mesh points for a given layoutType

        OUPUTS:
        -------
        The output is a list of two dictionaries as follows:
        [internalPoints,boundaryPoints], where the first contains
        the properties of the interior mesh points and the second that of the
        boundary points. The keys for each of these dictionaries are tuples
        containing the labels of the mesh points (a label is just a numbering of
        the mesh points) and the values are instasnces of hte class meshPoint. 
        """ 

        # Dictionaries whose keys are labels of the points in a 2-D grid and values
        # are an instance of the class meshPoint holding the informaiton about 
        # that mesh point
        self.boundaryPoints = {}
        self.internalPoints = {}

        # Rectangle 
        if self.layoutType.lower() == 'rectangle': 
            # Define the mesh for a rectanglular layout
            self.defineRectangleLayout()
        # Circle 
        elif self.layoutType.lower() == 'circle':
            # Define the mesh for a circular layout
            self.defineCircleLayout()

        return [self.internalPoints,self.boundaryPoints]


#--------- Sample implementation ------
if __name__ == "__main__":

    sq = layout2D(layoutType = 'rectangle', Lx = 1, Ly = 1.5, h = 0.2)
    [internalPoints,boundaryPoints] = sq.run()
    print '------- internalPoints ------------'
    for pointLabel in sorted(internalPoints.keys(),key=lambda tup:tup[0]):
        print 'label= ',pointLabel,'  (x,y) = (',internalPoints[pointLabel].x,',',internalPoints[pointLabel].y,')','  fW = ',internalPoints[pointLabel].fW,'  fE = ',internalPoints[pointLabel].fE,'  fN = ',internalPoints[pointLabel].fN,'  fS = ',internalPoints[pointLabel].fS,'  cW = ',internalPoints[pointLabel].coeffS,'  cE = ',internalPoints[pointLabel].coeffE,'  cN',internalPoints[pointLabel].coeffN,'  cS = ',internalPoints[pointLabel].coeffS,'\n'

    print '\n------- boundaryPoints ------------'
    for pointLabel in sorted(boundaryPoints.keys(),key=lambda tup:tup[0]):
        print 'label= ',pointLabel,'  (x,y) = (',boundaryPoints[pointLabel].x,',',boundaryPoints[pointLabel].y,')'

