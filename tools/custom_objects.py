import sys
from copy import copy, deepcopy
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

class customList(list):
    """
    This subclasses the python's built-in list. It does not allow the direct modification of
    the entries of the list (e.g., d[0] = 5 or del d[0]) and instead requires to use a 
    set method
    """
    def __setitem__(self,name, value):
        raise Exception('Direct modification of the entries of this list is not allowed. Use set_[property_name] method instead')

    def __delitem__(self,value):
        raise Exception('Direct modification of the entries of this list is not allowed. Use set_[property_name] method instead')

    def pop(self,value):
        raise Exception('Direct modification of the entries of this list is not allowed. Use set_[property_name] method instead')

    def append(self,value):
        raise Exception('Direct modification of the entries of this list is not allowed. Use set_[property_name] method instead')

    def __deepcopy__(self, memo):
        """
        Redfines the deeopcopy (original deepcopy does not work for this customized list class
        Source: http://stackoverflow.com/questions/37015123/deep-copying-a-user-defined-python-dictionary
        """
        def _deepcopy_dict(x, memo):
            y = {}
            memo[id(x)] = y
            for key, value in x.iteritems():
                y[deepcopy(key, memo)] = deepcopy(value, memo)
            return y
        cls = self.__class__
        result = cls.__new__(cls)
        result.__init__(_deepcopy_dict(self, memo))
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        return result

class customDict(dict):
    """
    This is a subclass python's built-in dictionary. It does not allow the direct modification of
    the entries of the dicitionary (e.g., d['a'] = 5 or del d['a']) and instead requires to use a 
    set method
    """
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __setitem__(self,key,value):
        raise Exception('Direct modification of the entries of this dictionary is not allowed. Use set_[property_name] method instead')

    def __delitem__(self,key):
        raise Exception('Direct modification of the entries of this dictionary is not allowed. Use set_[property_name] method instead')

    def __deepcopy__(self, memo):
        """
        Redfines the deeopcopy (original deepcopy does not work for this customized dict class
        Source: http://stackoverflow.com/questions/37015123/deep-copying-a-user-defined-python-dictionary
        """
        def _deepcopy_dict(x, memo):
            y = {}
            memo[id(x)] = y
            for key, value in x.iteritems():
                y[deepcopy(key, memo)] = deepcopy(value, memo)
            return y
        cls = self.__class__
        result = cls.__new__(cls)
        result.__init__(_deepcopy_dict(self, memo))
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        return result

#-------------------
class foo(object):
    def __init__(self,x, y):
        self.x = x
        self.y = y
    def __setattr__(self,attr_name,attr_value):
       print 'hello 0'
       if attr_name not in ['x','y']:
           self.__dict__[attr_name] = attr_value 
       if attr_name == 'x':
           print 'hello 1,  x = ',attr_value
           self.set_x(x = attr_value)
       if attr_name == 'y':
          self.set_y(y = attr_value)
    def set_x(self,x):
       print 'hello 2,  x = ',x
       self.__dict__['x'] = customList(x)
    def set_y(self,y):
       self.__dict__['y'] = customDict(y)
        
if __name__ == '__main__':
    f = foo(x = [1,2,3], y = {'a':1,'b':2,'c':[3,4]})
    print 'f.x = {} ,  f.y = {}'.format(f.x, f.y)
    print "f.x[1] = {} , f.y['a'] = {}".format(f.x[1],f.y['a'])
    #del f.x[1]
    #f.x[1] = 5
    #f.x.append(4)
    #f.x.pop(0)
    #f.x = [4,5]
    f.x += [4,5]
    print f.x

    #f.y['a'] = 5
    #f.y['d'] = 4
    #del f.x['a']
    #f.y = {'a':2, 'b':3}
    #print f.y

    z = deepcopy(f.y)
    z['c'].append(5) 
    print 'f.y = ',f.y
    print 'f.z',f.z

