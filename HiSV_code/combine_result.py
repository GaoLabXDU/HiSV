# read result file
from itertools import groupby
from collections import defaultdict


def combine(lst):
    pos = (j - i for i, j in enumerate(lst))
    t = 0
    for i, els in groupby(pos):
        l = len(list(els))
        el = lst[t]
        t += l
        yield range(el, el + l)


def group_1(data):
    last = data[0]
    start = end = 0
    for n in data[1:]:
        if n - last == 1:  # Part of the group, bump the end
            last = n
            end += 1
        else:  # Not part of the group, yield current group and start a new
            yield range(data[start], data[end] + 1)
            last = n
            start = end = end + 1
    # yield start, end
    yield range(data[start], data[end] + 1)


def group(chrlist):
    last = chrlist[0]
    start = end = 0
    for n in chrlist[1:]:
        if n == last:  # Part of the group, bump the end
            last = n
            end += 1
        else:  # Not part of the group, yield current group and start a new
            yield start, end
            last = n
            start = end = end + 1
    yield start, end


def group_position(pos1list, pos2list, binsize):
    result = defaultdict(list)
    l1 = sorted(set(pos1list), key=pos1list.index)
    pos1group = list(group_1(l1))
    for i in range(len(pos1group)):
        start_pos = pos1list.index(pos1group[i][0])
        end_pos = len(pos1list) - 1 - pos1list[::-1].index(pos1group[i][-1])
        cur_pos2_list = pos2list[start_pos:end_pos + 1]
        cur_pos1_list = pos1list[start_pos:end_pos + 1]
        cur_pos2_list_sort = sorted(cur_pos2_list)
        l2 = sorted(set(cur_pos2_list_sort), key=cur_pos2_list_sort.index)
        pos2group = list(group_1(l2))
        if len(pos2group) > 1:
            for j in range(len(pos2group)):
                new_pos1 = []
                for k in range(len(pos2group[j])):
                    if pos2group[j][k] in cur_pos2_list:
                        pos = cur_pos2_list.index(pos2group[j][k])
                        new_pos1.append(cur_pos1_list[pos])
                '''
                for m in range(len(pos1list)):
                    if pos2list[m] in pos2group[j]:
                        new_pos1.append(pos1list[m])

                for m in range(len(cur_pos1_list)):
                    if pos2list[m] in cur_pos2_list:
                        new_pos1.append(cur_pos1_list[m])
                '''
                new_pos1 = sorted(new_pos1)
                result['pos1_start'].append(new_pos1[0] * binsize)
                result['pos1_end'].append(new_pos1[-1] * binsize + binsize)
                result['pos2_start'].append(pos2group[j][0] * binsize)
                result['pos2_end'].append(pos2group[j][-1] * binsize + binsize)

        else:
            result['pos1_start'].append(pos1group[i][0] * binsize)
            result['pos1_end'].append(pos1group[i][-1] * binsize + binsize)
            result['pos2_start'].append(pos2group[0][0] * binsize)
            result['pos2_end'].append(pos2group[0][-1] * binsize + binsize)
    return result

'''
if __name__ == '__main__':

    filename = "/media/li/Data/T47D/hic_breakfinder_result/different_ratio_result/inter/T47D_ratio_20_result.txt"

    chr1, chr2, pos1, pos2 = read_filter_file(filename)
    group_chr = list(group(chr1))
    inter_result_file = "/media/li/Data/T47D/hic_breakfinder_result/different_ratio_result/inter/T47D_ratio_20_combine_result.txt"
    for i in range(len(group_chr)):
        cur_chr2 = chr2[group_chr[i][0]:group_chr[i][1]+1]
        no_order_pos1 = pos1[group_chr[i][0]:group_chr[i][1]+1]
        no_order_pos2 = pos2[group_chr[i][0]:group_chr[i][1]+1]
        sort_id = sorted(range(len(cur_chr2)), key=lambda k: cur_chr2[k])
        chr2list = [cur_chr2[i] for i in sort_id]
        cur_pos1 = [no_order_pos1[i] for i in sort_id]
        cur_pos2 = [no_order_pos2[i] for i in sort_id]
        print(chr2list)
        group_chr2 = list(group(chr2list))

        for j in range(len(group_chr2)):
            pos1list = cur_pos1[group_chr2[j][0]:group_chr2[j][1] + 1]
            pos2list = cur_pos2[group_chr2[j][0]:group_chr2[j][1] + 1]
            result = group_position(pos1list, pos2list)

            with open(inter_result_file, 'a') as out1:
                for k in range(len(result['pos1_start'])):
                    out1.write(str(chr1[group_chr[i][0]]) + '\t' + str(result['pos1_start'][k]) + '\t' + str(result['pos1_end'][k])
                              + '\t' + str(chr2list[group_chr2[j][0]]) + '\t' + str(result['pos2_start'][k]) + '\t' + str(result['pos2_end'][k]) + '\n')

'''

