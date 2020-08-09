
def get_sourcename(s1):
    s1name2 = get_dict_value(s1, 'name2', '')
    s1name = get_dict_value(s1, 'name', '')
    if s1name2 != '':
        s1name = s1name2
    return s1name

def mkdict(keys_array, values_array):
    m = {}
    for i, k1 in enumerate(keys_array):
        m[k1] = values_array[i]
    return m

def get_dict_value(dict1, key1, default):
    rst = default
    if key1 in dict1.keys():
        rst = dict1[key1]
    return rst

def convto_str_array(arr1):
    rst = []
    for a1 in arr1:
        rst.append(str(a1))
    return rst

def is_available(field):
    flag = True
    if get_dict_value(field, "is_available", True) is False:
        flag = False
    return flag


def merge_dict(dict1, dict2):
    return(dict2.update(dict1))


def conv_dict_from_arr(arr, header):
    print(len(arr), arr)
    print(len(header), header)
    d = {}
    for i in range(len(arr)):
        d[header[i]] = arr[i]
    return d
