import re


def get_bag_data(file_location):
    bags = {}
    with open(file_location, "r") as f:
        outer_bag_pattern = r"(^\w+ \w+ \w+)"
        inner_bag_pattern = r"(\d \w+ \w+ \w+)"
        for line in f.readlines():
            inner_bags = {}
            for inner_bag in re.findall(inner_bag_pattern, line):
                #print(inner_bag)
                bag_count = int(inner_bag[0])
                bag_name = inner_bag[2:] if bag_count == 1 else inner_bag[2:-1]
                inner_bags[bag_name] = bag_count
                #print(inner_bags)
            outer_bag_str = re.match(outer_bag_pattern, line).group(1)
            if outer_bag_str[-1] == "s":
                outer_bag = outer_bag_str[:-1]
            else:
                outer_bag = outer_bag_str
            bags[outer_bag] = inner_bags
        #print(bags)
        return bags


def part_a(file_location):
    bags = get_bag_data(file_location)

    holds_shiny_gold = {}

    def bag_search(bag_type):
        current_bag = bags[bag_type]
        if 'shiny gold bag' in current_bag:
            holds_shiny_gold[bag_type] = True
            return True

        current_bag_search = []
        for bag in current_bag:
            if holds_shiny_gold.get(bag):
                current_bag_search.append(holds_shiny_gold.get(bag))
            else:
                current_bag_search.append(bag_search(bag))
        holds_gold = any(current_bag_search)
        holds_shiny_gold[bag_type] = holds_gold
        return holds_gold

    for bag in bags:
        bag_search(bag)

    return sum(holds_shiny_gold.values())


def part_b(file_location):
    bags = get_bag_data(file_location)

    bag_count = {}

    def count_bags(bag_to_count):
        current_bag = bags[bag_to_count]
        inner_bag_count: list[int] = []
        for bag in current_bag:
            inner_bag_count.append(current_bag[bag])
            if bag_count.get(bag):
                inner_bag_count.append(bag_count[bag] * current_bag[bag])
            else:
                inner_bag_count.append(count_bags(bag) * current_bag[bag])
        total_bags = sum(inner_bag_count)
        bag_count[bag_to_count] = total_bags
        return total_bags

    shiny_gold_bag_count = count_bags('shiny gold bag')
    return shiny_gold_bag_count

get_bag_data('advantcode_day7_data.txt')

#a = part_a('advantcode_day7_data.txt')
#print(a)
#print(part_b('advantcode_day7_data.txt'))