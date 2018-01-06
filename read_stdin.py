import sys

"""To run the script, enter the following command:
    $ read_stdin.py < example_input.txt
"""

if __name__ == '__main__':
    lines = sys.stdin.read().splitlines()
    for index, line in enumerate(lines):
        print(f'Line {index} is {len(line)} characters long.')

    