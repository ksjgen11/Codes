# Day 6_part 2

def read_group_files(filename):
    '''read strings from file'''
    file_open = open(filename, 'r')
    answers = file_open.read().splitlines()

    keys = []
    item = []
    items = []
    total_answers = {}
    x = 0
    for answer in answers:
        if answer == '':
            keys.append(x)
            x += 1
            if item:
                items.append(item)
                item = []
        else:
            item.append(answer)
    keys.append(x)
    items.append(item)  # last sequence appending

    for i in range(len(keys)):
        total_answers[keys[i]] = sorted(items[i])

    return total_answers

answers = read_group_files('advantcode_day6_data.txt')
#print(answers)

def find_overlap_number(dictionary_input):
    answers = dictionary_input
    count = {}
    number = 0
    for keys, values in answers.items():
        count[keys] = {}
        for elements in values:
            for element in elements:
                try:
                    count[keys][element] += 1
                except:
                    count[keys][element] = 1
        #print(count)
        for dic in count[keys]:
            if count[keys][dic] == len(answers[keys]):
                number += 1
                #print(count[keys][dic], len(answers[keys])

    return number



print(find_overlap_number(answers))


