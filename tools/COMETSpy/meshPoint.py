from __future__ import division
import math
sys.path.append('../../')
from tools.customError import *


class meshPoint(object):
    """
    This class holds the properties of a mesh point
    """

    def __init__(self,type,x,y,z = None ,xLabel = None, yLabel = None,zLabel = None, fE = None, fW = None, fN = None, fS = None,metabs = None, species = None):

        """
        INPUTS:
        -------
                         type: Type of the mesh point. Eligible inputs are 'boundary' 
                               and 'interior'
                        x,y,z: x-, y- and z-coordinates of the mesh point (z-coordinate 
                               is optional)
        xlable, ylabel,zlabel: x-, y- and z-label of the mesh point (z-label is optional)
                               labels can be anything and are chosen by tthe user. 
                               We typically number the points and assign these numbers to
                               them as lables
             fE,fW, fN and fS: Distance from the points located on the east, west, south 
                               and north divided by the mesh size (h)
                       metabs: A dictionary whose names ar the names of the metabolites
                               presents in this mesh point and values are instances of the
                               class sharedMetab 
                    organisms: A dictionary whose names ar the names of the microorganisms
                               presents in this mesh point and values are their concentrations
        """

        # Type of the mesh point (boundary or interior)
        if type.lower() not in ['boundary','internal']:
            raise customError("**Error! Invalid mesh point type (allowable choices for the mesh point 'type' are 'boundary' and 'interior')")
        else:
            self.type = type

        # x-coordinate
        self.x = x

        # y-coordinate
        self.y = y

        # z-coordinate
        self.z = z

        # x-axis component of the mesh point label
        self.xLabel = xLabel

        # y-axis component of the mesh point label
        self.yLabel = yLabel

        # y-axis component of the mesh point label
        self.zLabel = zLabel

        # Distance from the points located on the east, west, south and north
        # divided by self.h
        self.fE = fE
        self.fW = fW
        self.fN = fN
        self.fS = fS

        # Calculate the coefficients needed to compute the Laplacian if the 
        # values of fE, fW, fN and fS are provided
        if self.fE is not None and self.fW is not None and self.fN is not None and self.fS is not None:
            self.LapCoeffCalc()

    def LapCoeffCalc(self):
        """
        A function to calculate the coefficients needed to compute the 
        Laplacian (d^C/dx^2, d^2C/dy^2, ...) 
        """
        # East
        self.cE = 2/(self.fE*(self.fW + self.fE)) 

        # West
        self.cW = 2/(self.fW*(self.fW + self.fE)) 

        # North
        self.cN = 2/(self.fN*(self.fS + self.fN)) 

        # South
        self.cS = 2/(self.fS*(self.fS + self.fN)) 

        # cX
        self.cX = 1/(self.fE*self.fW)

        # cY
        self.cY = 1/(self.fS*self.fN)


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

