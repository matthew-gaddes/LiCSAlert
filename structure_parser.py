#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:37:28 2024

@author: matthew
"""

print("started")

import os
import ast
import matplotlib.pyplot as plt
import sys

import pdb

from pathlib import Path

def read_py_file_lines(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines



#%% 
class py_function:
    
    def __init__(self, name, file, line_start):
        self.name = name
        self.file = file
        self.line_start = line_start
        
    def find_line_end(self):
        """ Find the line that the function ends on.  
        """
        
        lines = read_py_file_lines(self.file)
        
        end_found = False
        line_n = self.line_start + 1
        while (not end_found):
            line = lines[line_n]
            # function ends either when the next function starts
            if line[:3] == "def":
                self.line_end = (line_n -1) 
                end_found = True
            
            # or the file ends, for the last function.  
            if line_n == (len(lines) -1):
                self.line_end = (line_n -1) 
                end_found = True
            line_n += 1
            
                
#%%


class TreeNode:
    """ To store a heirarchical set of py_function to describe 
    all the function calls within a Python package.  """
    
    def __init__(self, value):
        self.value = value
        self.children = []

    def add_child(self, child_node):
        self.children.append(child_node)

    def __repr__(self, level=0):
        """ Returns a py_function.name attribute to be readable.  
        """
        # each level of function is indented by one tab (\t)
        ret = "\t" * level + repr(self.value.name) + "\n"
        for child in self.children:
            ret += child.__repr__(level + 1)
        return ret  
    
    def find_function_n(function_name, py_functions):
        """
        """
        for func_n, py_f in enumerate(py_functions):
            if py_f.name == function_name:
                return func_n
        raise Exception(f"Failed to find {function_name} in py_functions")
        
        
#%% 


def find_all_functions(directory, verbose = True):
    """ Find all the functions in a python directory (including child 
    directories)
    """

    # Find all the functions in a directory of python files.  
    py_functions = []
    # recursively look at all directories
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".py"):
                print(f"Opening {file}")
                file_path = Path(root) / file
                lines = read_py_file_lines(file_path)

                for line_n, line in enumerate(lines):
                    # if this line is the start of a function
                    if line[:3] == "def":
                        # get the name (between def and first ( )
                        function_name = line.split('(')[0][4:]
                        # add to list of functions
                        py_functions.append(py_function(function_name, file_path,
                                                            line_n))
                # debug
                # for py_f in py_functions:
                #     print(py_f.name)
    
    # get the line end numbers for the functions
    for py_f in py_functions:
        py_f.find_line_end()
        
    # single list of all the function names
    py_function_names = [py_f.name for py_f in py_functions]
    
    if verbose:
        for py_f in py_functions:
            print(f"{py_f.name} : {py_f.line_start} - {py_f.line_end} ")
        
    return py_functions, py_function_names

#%%

def find_function_calls(script, line_start, line_end,
                        functions_of_interest):
    """ Find the functions called by a function
    """
    import ast
    
    with open(script, "r") as file:
        tree = ast.parse(file.read(), filename=script)


    # 1: get all the functions used in that block of code (line_start - line_end)    
    all_function_calls = []
    class FunctionCallVisitor(ast.NodeVisitor):
        def visit_Call(self, node):
            if isinstance(node.func, ast.Name):
                if line_start <= node.lineno <= line_end:
                    all_function_calls.append(node.func.id)
            self.generic_visit(node)
    
    FunctionCallVisitor().visit(tree)
    
    # 2: reject any that aren't the ones we're interested in
    # i.e. usually things that aren't in our module/package
    
    function_calls = []
    for function_call in all_function_calls:
        if function_call in functions_of_interest:
            function_calls.append(function_call)
    return function_calls

#%%

def find_function_n(function_name, py_functions):
        """
        """
        for func_n, py_f in enumerate(py_functions):
            if py_f.name == function_name:
                return func_n
        raise Exception(f"Failed to find {function_name} in py_functions")


#%%


def add_functions_to_tree(func_tree_node, py_functions, py_function_names):
    """
    """

    # find all functions used by the function.  
    func_calls = find_function_calls(func_tree_node.value.file,
                                     func_tree_node.value.line_start, 
                                     func_tree_node.value.line_end, 
                                     py_function_names)
    
    # if there are no functions of interest in this function, nothing more:
    if len(func_calls) == 0:
        return 

    # else add those functions as children to this node.  
    else:
        print(f"In {func_tree_node.value.name} the following functions were found: ")
        for func_call in func_calls:
            print(f"    {func_call}")
        
        # convert the strings back to py_functions, and add as children to node.          
        for child_n, func_call in enumerate(func_calls):
            index_n = find_function_n(func_call, py_functions)    
            py_func = py_functions[index_n]
            func_tree_node.add_child(TreeNode(py_func))
            
            # Recursive look for child functions of this function.  
            print(f"Looking for functions that are called in {py_func.name}")
            add_functions_to_tree(func_tree_node.children[child_n], py_functions,
                                  py_function_names)

        return func_tree_node


#%% 

def build_function_call_tree(directory, root_function):
    """ Build a heirarchical tree of all the functions in a python
    package.  
    
    """

    # find all the functions in the directory (names, files, line numbers)
    py_functions, py_function_names = find_all_functions(directory, verbose = True)
          
    # Get the root function (i.e. the one called)
    func_n = find_function_n(root_function, py_functions)    
    root_function = py_functions[func_n]

    # start the heirarchical tree
    func_tree_node = TreeNode(root_function)

    # populate the heirarchical tree recursively.  
    func_tree = add_functions_to_tree(func_tree_node, py_functions, py_function_names)
        
    return func_tree


#%%    

directory = "/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/00_LiCSAlert_GitHub"
root_function = "LiCSAlert_monitoring_mode"

func_tree = build_function_call_tree(directory, root_function)


func_tree

#%%




sys.exit()




def print_tree(tree, indent=0):
    for key, value in tree.items():
        print("  " * indent + key)
        print_tree(value, indent + 1)













def plot_tree(tree, root_function):
    fig, ax = plt.subplots(figsize=(8, 12))  # Make the figure taller and narrower
    ax.axis('off')

    def draw_tree(tree, x=0.5, y=1, depth=1):
        if not tree:
            return x, y
        
        node = list(tree.keys())[0]
        ax.text(x, y, node, ha='center', va='center', bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="skyblue"))

        children = list(tree[node].keys())
        if not children:
            return x, y
        
        new_y = y - 1 / depth
        new_x = x - (len(children) - 1) / (2 * depth)
        for i, child in enumerate(children):
            ax.plot([x, new_x + i / depth], [y, new_y], 'k-', lw=1)
            draw_tree({child: tree[node][child]}, new_x + i / depth, new_y, depth * 2)
        
        return x, y

    draw_tree(tree, x=0.5, y=1, depth=1)
    plt.title('Function Call Tree')
    plt.show()

# Usage example:
directory_path = Path("/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/00_LiCSAlert_GitHub/licsalert")
root_function_name = "LiCSAlert_monitoring_mode"

function_tree = build_function_call_tree(directory_path, root_function_name)
print_tree(function_tree)
plot_tree(function_tree, root_function_name)
