
def is_available(field):
    flag = True
    if 'is_available' in field.keys() and field['is_available'] == False:
        flag = False
    return flag


