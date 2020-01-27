class Node:
    def __add__(self, other):
        return Sum(self, other)


class Sum(Node):
    def __init__(self, left, right):
        self._left = left 
        self._right = right 

    @property
    def children(self):
        return (self._left, self._right)

class Constant(Node):
    def __init__(self, value):
        self._value = value 

    children = ()


import random 

def tree(N):
    root = Constant(1) + Constant(1)
    for i in range(N):
        if random.random() < 0.5:
            root = root + Constant(1)
        else:
            root = Constant(1) + root 
    return root 


def count_nodes(tree):
    stack = [tree] 
    total = 0

    while len(stack) > 0:
        e = stack.pop()
        total += 1
        for child in e.children:
            stack.append(child)

    return total

print(count_nodes(tree(100)))