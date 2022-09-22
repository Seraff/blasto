#!/usr/bin/env python3

from numpy import full


def get_anticodon_real_seq(full_string, start, end):
    as_list = list(full_string)
    as_list.insert(end+1, ']')
    as_list.insert(start, '[')
    print(''.join(as_list))

    before = full_string[start-2:start]
    after = full_string[end+1:end+3]

    if before == '>.':
        start -= 2

        while full_string[start] == '>':
            start -= 1
        start += 1

    if after == '.>':
        end += 2

        while full_string[end] == '>':
            end += 1
        end -= 1

    return full_string[start:end+1]


seq = '..>>>.>..>>>>..'
start = 6
end = 6

new_str = get_anticodon_real_seq(seq, start, end)
print(new_str)
