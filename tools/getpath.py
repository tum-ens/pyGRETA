import os
path = os.path.abspath(__file__)
print(path)
print(os.getcwd())

print(os.listdir('..'))

home = os.path.expanduser("~")
print(os.listdir(home + "/pygreta"))