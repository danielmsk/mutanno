
def get_dict_value(dict1, key1, default):
    rst = default
    if key1 in dict1.keys():
        rst = dict1[key1]
    return rst


def is_available(field):
    flag = True
    if get_dict_value(field, "is_available", True) == False:
        flag = False
    return flag

def merge_dict(dict1, dict2):
	return(dict2.update(dict1)) 