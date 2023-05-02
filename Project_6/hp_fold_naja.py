

hp = "hhppppphhppphppphp"

def even_odd(string): 
    list_index = []
    for i in range(len(string)): 
        if string[i] == "h":
            list_index.append(i)
    return list_index 

print(even_odd(hp))



def score_energy_left(even_odd_list):
    energy_score = 0
    energy_list = []
    low = 0 
    up = 1
    
    while even_odd_list[low] < even_odd_list[-up] -1:
        if even_odd_list[low] % 2 == 0: 
            if even_odd_list[-up] % 2 == 1: 
                energy_score += 1
                energy_list.append(even_odd_list[low])
                energy_list.append(even_odd_list[-up])
                up += 1
                low += 1
            else: 
                up += 1

        elif even_odd_list[low] % 2 == 1: 
            if even_odd_list[-up] % 2 == 0: 
                energy_score += 1
                energy_list.append(even_odd_list[low])
                energy_list.append(even_odd_list[-up])
                up += 1
                low += 1
            else: 
                up += 1

    print(energy_list, "e")    
    return energy_score
                    


def score_energy_right(even_odd_list):
    energy_score = 0
    energy_list = []
    low = 0 
    up = 1
    while even_odd_list[low] < even_odd_list[-up] -1:
        if even_odd_list[-up] % 2 == 0: 
            if even_odd_list[low] % 2 == 1: 
                energy_score += 1
                energy_list.append(even_odd_list[low])
                energy_list.append(even_odd_list[-up])
                up += 1
                low += 1
            else: 
                low += 1

        if even_odd_list[-up] % 2 == 1: 
            if even_odd_list[low] % 2 == 0: 
                energy_score += 1
                energy_list.append(even_odd_list[low])
                energy_list.append(even_odd_list[-up])
                up += 1
                low += 1
            else: 
                low += 1
    print(energy_list, "e")          
    return energy_score
                    
print(score_energy_left(even_odd(hp)), 1)
print(score_energy_right(even_odd(hp)), 1)

hp_2 = "hphphhhppphhhhpphh"
print(even_odd(hp_2))
print(score_energy_left(even_odd(hp_2)), 2)
print(score_energy_right(even_odd(hp_2)), 2)


hp_3 = "phpphphhhphhphhhhh"
print(even_odd(hp_3))
print(score_energy_left(even_odd(hp_3)), 3)
print(score_energy_right(even_odd(hp_3)), 3)


hp_4 = "hphpphhphpphphhpphph"
print(even_odd(hp_4))
print(score_energy_left(even_odd(hp_4)), 4)
print(score_energy_right(even_odd(hp_4)), 4)