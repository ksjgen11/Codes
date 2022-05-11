
'''
Different ways of concatenating strings using loops and joins
BTW: Note that multi-line comments are done using three 's before and after
'''

# Function for Concatenation of numbers using a loop 
def loop_concat(numbers):
    string = ''
    for p in numbers:
        string += str(p)
    return string

# Function for Concatenation of numbers using Functional style 
def join_concat(numbers):
    string = ''.join(map(str,numbers))
    return string
