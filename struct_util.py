
def get_dict_value(dict1, key1, default):
    rst = default
    if key1 in dict1.keys():
        rst = dict1[key1]
    return rst


def is_available(field):
    flag = True
    if 'is_available' in field.keys() and field['is_available'] == False:
        flag = False
    return flag


