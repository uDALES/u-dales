import os
import netCDF4
import numpy as np
import pickle


class Exp(object):
    """Class for uDALES simulation objects.
    By default this class has 3 groups to 
    store content: parameters, data, stats.
    Initialise object with 'Exp(name)'."""

    def __init__(self, name):
        self.name = name
        self.parameters = {}
        self.data = {}
        self.stats = {}


    def groups(self):
        groups = []
        for key, val in vars(self).items():
            if isinstance(val, dict):
                groups.append(val)
        return groups


    def __getattr__(self, name):
        for group in self.groups():
            if name in group:
                return group[name]
        return self.__dict__[name]


    def __setattr__(self, name, value):
        for group in self.groups():
            if name in group:
                group[name] = value
                return
        self.__dict__[name] = value


    def add(self, name, value, group=None):
        """Add a new entry to the class.
        A group can be optionally specified."""
        if group is None:
            self.__dict__[name] = value
        elif group.casefold() in self.__dict__:
            self.__dict__[group.casefold()][name] = value
        else:
            print('Could not add item.')

    
    def delete(self, name, group=None):
        """Delete an entry from the class.
        If the entry is part of a group, then 
        the group needs to be specified."""
        if group is None:
            del self.__dict__[name]
        elif group.casefold() in self.__dict__ and name in self.__dict__[group.casefold()]:
            del self.__dict__[group.casefold()][name]
        else:
            print('Could not delete item.')


    def save(self, path):
        """Save object to file."""
        # create directory if it does not exist
        directory = os.path.dirname(path)
        os.makedirs(directory, exist_ok=True)
        with open(path, 'wb') as file:
            ExpStored.save(self, file)


    def show(self, level='content', group=None):
        """List all the object's content.
        The level of depth can be specified with the
        'level' keyword: 'group', 'keys', 'content'.
        To only show the content from one group, use
        the 'group' keyword."""
        if level.lower() in ['group', 'groups', 'g', 1]:
            depth = 1
        elif level.lower() in ['key', 'keys', 'k', 2]:
            depth = 2
        elif level.lower() in ['content', 'all', 'c', 'a', 3]:
            depth = 3
        else:
            raise KeyError('Level not found. Choose from: "group", "keys", "content".')

        if group is None:
            content = self.__dict__
        elif group.casefold() in self.__dict__:
            content = self.__dict__[group.casefold()]
            print(group.upper())
        else:
            raise KeyError('Could not find group. Try group=None.')

        for key, val in content.items():          
            if isinstance(val, dict):
                if depth == 1:
                    print(key)
                elif depth > 1:
                    print(key.upper())

                    for subkey, subval in val.items():
                        if depth == 2:
                            print("\t", subkey)
                        elif depth == 3:
                            if isinstance(subval, (str, int, float, complex)):
                                print("\t", subkey, ":", subval)
                            else:
                                print("\t", subkey, ":", type(subval))

            else:
                if depth == 2:
                    print(key)
                elif depth == 3:
                    if isinstance(val, (str, int, float, complex)):
                        print("\t", key, ":", val)
                    else:
                        print("\t", key, ":", type(val))

    
class ExpStored(object):
    """Class for saving uDALES simulation objects 
    as pickle files."""

    def __init__(self):
        pass
    
    @staticmethod
    def create_type(name='DynamicType', dict={}):
        return type(name, (object,), dict)
    
    @staticmethod
    def save(t, fh):
        dict = t.__dict__.copy()
        name = t.name
        for key in dict.keys():
            if key.startswith('__') and key.endswith('__'):
                del dict[key]
        pickle.dump((name, dict), fh)
        
    @classmethod
    def load(cls, fh):
        name, dict = pickle.load(fh)
        return cls.create_type(name, dict)

    
def load_saved_exp(file):
    """Load a saved object file and create a 
    new object with its content."""
    with open(file, 'rb') as f:
        exp_saved = ExpStored.load(f) 
    exp = Exp(exp_saved.name)  # copy to new exp object
    for name, value in exp_saved.__dict__.items():
        if not name.startswith("__"):
            exp.__setattr__(name, value)
    return exp


def netcdf_handling(exp, maxdim=2):
    """Handles netcdf files, which cannot be saved to
    a pickle file. Converts the netcdf data to numpy 
    arrays for netcdf variables with dimensions less or 
    equal to 'maxdim'. Deletes netcdf variables from object
    for dimensions greater than 'maxdim'."""
    netcdfs_add = []
    netcdfs_delete = []
    for key, val in exp.data.items():
        if isinstance(val, netCDF4._netCDF4.Variable):
            if len(val.shape) > maxdim:
                netcdfs_delete.append(key)
            else:
                netcdfs_add.append(key)
    [exp.delete(d, 'data') for d in netcdfs_delete]
    [exp.add(a, exp.data[a][:], 'data') for a in netcdfs_add]
    print("Added as numpy arrays:\n", netcdfs_add)
    print("Removed from exp object:\n", netcdfs_delete)
    return exp
