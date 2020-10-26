import skrf as rf
from matplotlib import pyplot as plt
import numpy as np
import math as math
import statistics
import random
import toolz
import ttg
import itertools
import pickle
import pathlib
import os
#returns the names of the files in the directory data as a list


###### Подготовильтельный блок
list_of_files = os.listdir("dsa_almaz_s2p_file"); # директория

CellPath = os.listdir("dsa_almaz_s2p_file");

Num_Of_Bits=int(len(list_of_files)/2); # Определили число битов в комбинации;
#Создаем таблицу истинности
# Лист битов ['Bit0', 'Bit1', 'Bit2', 'Bit3', 'Bit4', 'Bit5']
BitList = [];
for i in range(0,Num_Of_Bits,1):
    BitList.append(str('Bit'+str(i)));
#Таблица истинности
Bit = ttg.Truths(BitList, ascending=True) # развенутная , младший  - бит последний
###### Таблица Истиности готова;
print(Bit);
number_of_states = 2**Num_Of_Bits;  ## число состояний
#Надо разобрать еще все файлы, 
Cell_Name_List = [];
Cell_P1dB_List = [];
Cell_Truth_List = [];

################################
PhaseStep = 360/(2**Num_Of_Bits);
NormalPhaseList= [];
for state in range(number_of_states):
    Phase = state*PhaseStep;
    NormalPhaseList.append(Phase);
#####################################
NormasAttenuationList = [];
AttenuationStep = -0.5;
for state in range(number_of_states):
    Attenuation = state*AttenuationStep;
    NormasAttenuationList.append(Attenuation);
################################################
################################################
for i in range(len(list_of_files)):
    Split_String =  list_of_files[i].split("_");
    Cell_Name_List.append(Split_String[0]);
    Cell_P1dB_List.append(Split_String[1]);
    Cell_Truth_List.append(Split_String[2][:1]);

print(BitList)
print(BitList)




class Cell:  # классс ячеек - 2* число бит
    def __init__(self, CellValue, P1dB, ON_OFF, Network, Cell_S21_MF_dB):
            self.CellValue = CellValue;
            self.P1dB = P1dB;              
            self.ON_OFF = ON_OFF; 
            self.Network = Network;###########################
            self.Cell_S21_MF_dB = Cell_S21_MF_dB;
   
def GetCells(CellPath):  # Создание всех ячееек
 
    Cell_Name_List = [];
    Cell_P1dB_List = [];
    Cell_Truth_List = [];
    list_of_files = CellPath;

    for i in range(len(list_of_files)):
        Split_String =  list_of_files[i].split("_");
        Cell_Name_List.append(Split_String[0]);
        Cell_P1dB_List.append(Split_String[1]);
        Cell_Truth_List.append(Split_String[2][:1]);
    CellList = [];

    for i in range(len(list_of_files)):
        CellValue = Cell_Name_List[i];
        P1dB = Cell_P1dB_List[i];
        ON_OFF = Cell_Truth_List[i];
        Network = rf.Network('dsa_almaz_s2p_file/'+str(CellValue)+'_'+str(P1dB)+'_'+str(ON_OFF)+'.S2P');
        #Network = rf.Network('to Vuchich/'+str(CellValue)+'_'+str(P1dB)+'_'+str(ON_OFF)+'.S2P');

        Frequency = Network.f.tolist(); # лист частот
        S21_dB = [];
        for f in range(len(Frequency)):    
            S21_dB.append(float(Network.s21[str(Frequency[f])+'hz'].s_db[...]));         
        MF = math.ceil(len(Frequency)/2);  # центральная ччастотта
        Cell_S21_MF_dB  = S21_dB[MF]; # выбор индекса на центральной частотте

        Temp_Cell = Cell(CellValue, P1dB, ON_OFF, Network, Cell_S21_MF_dB);
        CellList.append(Temp_Cell);
    return CellList;

#Лист всех ячеек
List_Cells = GetCells(CellPath);
############
UniqueNames = list(set(Cell_Name_List));
UniqueNames.sort(key=float); ## отсортироваали
UniqueNames.reverse(); ## развернули, младший бит последний
print(UniqueNames)
#########################################




def GetStateData(StateNumber):
    ##### Надо это обернуть в функцию
    CellList = []; # лист обьектов ячейки Сell для состояния StateString
    CascadeList =  []; # лист 2 полюсников rf.Network для состояния, в дальнейшем создания системы и перестановокж
    # Надо заполнить CascadeList по порядку в соотвествиие с таблицей истинност
    # Значит по порядку, это перебор по именам UniqueNames
    StateString =  StateNumber; # 000 000
    #StateString = 1; # 000 001 / младший бит последний
    #....
    #StateString = 63; # 111 111
    for i in range(len(UniqueNames)):
#    # теперь идем по списку ячеек и выбираем нужные
        for cell in range(len(List_Cells)):
            SelectCell = List_Cells[cell];
            if (SelectCell.CellValue==UniqueNames[i] and SelectCell.ON_OFF==str(int(Bit.base_conditions[StateString][i])) and SelectCell not in CellList):
                CellList.append(SelectCell)
                CascadeList.append(SelectCell.Network);
            else:
                continue;
##################### Return состояние - лист каскадов
    return CellList, CascadeList;


###########################   Проверка вроде работает
CellList, CascadeList= GetStateData(13);
for states in range(len(CellList)):
    print(CellList[states].CellValue+' '+CellList[states].ON_OFF);
##############################################

class State:
     def __init__(self, Number, CellList, CascadeList, 
                  P1dB,  # точка сжатия для данного состояния при каскадном включении
                  
                  S11_dB,  # лист значений
                  S22_dB,  # /-/-/-/-/-/-/-
                  S21_dB,  # /-/-/-/-/-/-/-
                  S21_deg, # /-/-/-/-/-/-/-

                  S11_dB_Max, # максимальное значение из списка частот
                  S22_dB_Max, # максимальное значение из списка частот

                  S21_dB_MF, # дБ на центральной частоте
                  S21_deg_MF, # Градус на центральной частоте
                  CombinationOrder): 
            
            self.Number = Number;
            self.CellList = CellList;
            self.CascadeList = CascadeList;
            self.P1dB = P1dB;

            self.S11_dB = S11_dB;
            self.S22_dB = S22_dB;
            self.S21_dB = S21_dB;
            self.S21_deg = S21_deg;

            self.S11_dB_Max = S11_dB_Max;
            self.S22_dB_Max = S22_dB_Max;

            self.S21_dB_MF = S21_dB_MF;
            self.S21_deg_MF = S21_deg_MF;
            self.CombinationOrder = CombinationOrder; #имя комбинации





def GetP1dBCascade(CellList):
    P1dB_Value_List = [];  # лист значений для данной комбинации
    S21_Value_List  = []; 

    for c in range(len(CellList)):
        P1dB_Value_List.append(CellList[c].P1dB);
        S21_Value_List.append(CellList[c].Cell_S21_MF_dB);
 
    #Cчитаем тотал P1dB
    #Переводим п1дб и s21 в разы;
    P1dB_Value_List_R = [];
    S21_Value_List_R = [];
    for i in range(len(P1dB_Value_List)):
        P1dB_Value_List_R.append(10**(float(P1dB_Value_List[i])/10));
        S21_Value_List_R.append(10**(float(S21_Value_List[i])/10));
    ###Дальше считаем
    #Делаем цикл по значением в листе P1dB_Value_List_R а для S21_Value_List_R -1
    SummList = [];
    GainFactArray = [];

    first = 1/float(P1dB_Value_List_R[0]); # Всегда
    GainMult = float(1); # временная переменная
    
    for x in range(len(S21_Value_List_R)-1):
        Gain = float(S21_Value_List_R[x]); #Cейчас
        GainMult = GainMult*Gain;       # следующая
        GainFactArray.append(GainMult);
        
    SuperSumma = 0;
    for cascade in range(1,len(P1dB_Value_List)):
        SuperSumma = SuperSumma + (GainFactArray[cascade-1]/P1dB_Value_List_R[cascade]);

    P1dBTotal = (SuperSumma+first)**(-1);
    P1dBTotal = 10*math.log10(P1dBTotal);
    P1dB = P1dBTotal;

    return P1dB;










def GetState(StateNumber, PermutaionNumber):
    CellList, CascadeList= GetStateData(StateNumber);  ## получили данные для состояния
    Cascade = CascadeList;
    Cascade_permutations = list(itertools.permutations(Cascade, len(BitList))); # 720 лист этих каскадов - 1 обьект листа это лист с 6 ячейками он же будет одинаково выдавать если состояние не менял??
    Cell_Permutaitions = list(itertools.permutations(CellList, len(BitList)))  # 720 состояний ячеек, отсюда можно взять последовательность
    CombinationOrder = [];
    for i in range(len(CellList)):
        CombinationOrder.append(Cell_Permutaitions[PermutaionNumber][i].CellValue);
    Number = StateNumber; # номер состояния
    Network = rf.cascade_list(Cascade_permutations[PermutaionNumber]);  # Берем все данные оттуда из собранного каска да
    Frequency = Network.f.tolist(); # лист частот
    S11_dB = [];   ### точки в по частотам
    S22_dB = [];
    S21_dB = [];
    S21_deg = [];
    for f in range(len(Frequency)):
        S11_dB.append(float(Network.s11[str(Frequency[f])+'hz'].s_db[...]));
        S21_dB.append(float(Network.s21[str(Frequency[f])+'hz'].s_db[...]));
        S22_dB.append(float(Network.s22[str(Frequency[f])+'hz'].s_db[...]));
        S21_deg.append(float(Network.s21[str(Frequency[f])+'hz'].s_deg[...]))
        #print('waiting');
    S11_dB_Max = max(S11_dB);
    S22_dB_Max = max(S22_dB);
    MF = math.ceil(len(Frequency)/2);  # центральная ччастотта
    S21_dB_MF  = S21_dB[MF]; # выбор индекса на центральной частотте
    S21_deg_MF = S21_deg[MF]
    P1dB = GetP1dBCascade(CellList);  # обернул в функцию
    return State(Number, CellList, CascadeList, 
                  P1dB,  # точка сжатия для данного состояния при каскадном включении       
                  S11_dB,  # лист значений
                  S22_dB,  # /-/-/-/-/-/-/-
                  S21_dB,  # /-/-/-/-/-/-/-
                  S21_deg, # /-/-/-/-/-/-/-
                  S11_dB_Max, # максимальное значение из списка частот
                  S22_dB_Max, # максимальное значение из списка частот
                  S21_dB_MF, # дБ на центральной частоте
                  S21_deg_MF, # Градус на центральной частоте)
                  CombinationOrder);


State_0 = GetState(0,0);

for cell in range(len(State_0.CellList)):
    print(str(State_0.CellList[cell].CellValue)+' '+str(State_0.CellList[cell].ON_OFF));

print('sho');


UniqueNames = UniqueNames;

def GetStateDataNew(StateNumber, UniqueNames):
    ##### Надо это обернуть в функцию
    CellList = []; # лист обьектов ячейки Сell для состояния StateString
    CascadeList =  []; # лист 2 полюсников rf.Network для состояния, в дальнейшем создания системы и перестановокж
    # Надо заполнить CascadeList по порядку в соотвествиие с таблицей истинност
    # Значит по порядку, это перебор по именам UniqueNames
    StateString =  StateNumber; # 000 000
    #StateString = 1; # 000 001 / младший бит последний
    #....
    #StateString = 63; # 111 111
    for i in range(len(UniqueNames)):
#    # теперь идем по списку ячеек и выбираем нужные
        for cell in range(len(List_Cells)):
            SelectCell = List_Cells[cell];
            if (SelectCell.CellValue==UniqueNames[i] and SelectCell.ON_OFF==str(int(Bit.base_conditions[StateString][i])) and SelectCell not in CellList):
                CellList.append(SelectCell)
                CascadeList.append(SelectCell.Network);
            else:
                continue;
##################### Return состояние - лист каскадов
    return CellList, CascadeList;



S = 1
print('****************************************************');
UniqueNames = UniqueNames;

CellList0, CascadeList0 = GetStateDataNew(S,UniqueNames)

for i in range(len(CellList0)):
    print(str(CellList0[i].CellValue)+'   '+str(CellList0[i].ON_OFF))
##############################################################
a = 0;
b = 5;
UniqueNames[a], UniqueNames[b] = UniqueNames[b], UniqueNames[a]
print('****************************************************');

CellList0, CascadeList0 = GetStateDataNew(S,UniqueNames)

for i in range(len(CellList0)):
    print(str(CellList0[i].CellValue)+'   '+str(CellList0[i].ON_OFF))






####Вроде все работает.
###################### Делаем систему состояний

#class StateSystemDevice:
#     def __init__(self, CombinationOrder,

#                        List_Of_States,
#                        List_Of_Max_S11_dB,
#                        List_Of_Max_S22_dB,
#                        List_Of_MF_P1dB,

#                        List_Of_MF_S21_dB,
#                        List_Of_MF_S21_deg,
                                                
#                        Max_S11_dB,
#                        Max_S22_dB,
#                        Min_P1dB,
                        
#                        RMS_S21_dB_DSA,
#                        RMS_S21_deg_DSA,

#                        RMS_S21_dB_PhSh,
#                        RMS_S21_deg_PhSh,
                        
#                        FitnessValue_DSA,
#                        FitnessValue_PhSh): 
            
#            self.CombinationOrder = CombinationOrder;   #Комбинация

#            self.List_Of_States = List_Of_States;   # лист непосредственно состояний
#            self.List_Of_Max_S11_dB = List_Of_Max_S11_dB;  # лист максимальных значений для всех состояний
#            self.List_Of_Max_S22_dB = List_Of_Max_S22_dB; # лист максимальных значений для всех состояний
#            self.List_Of_MF_P1dB = List_Of_MF_P1dB;      # лист П1дБ значений для всех состояний

#            self.List_Of_MF_S21_dB = List_Of_MF_S21_dB;   # Лист с21 дБ для расчета СКО
#            self.List_Of_MF_S21_deg = List_Of_MF_S21_deg;  # Лист с21 град для расчета СКО

#            self.Max_S11_dB = Max_S11_dB;                   #Выбор макисмального значения из листа List_Of_Max_S11_dB  - Для фитнесс функции
#            self.Max_S22_dB = Max_S22_dB;            #Выбор макисмального значения из листа List_Of_Max_S22_dB   - - Для фитнесс функции

#            self.Min_P1dB = Min_P1dB;               #Выбор минимального значения из листа List_Of_MF_P1dB;   -- Для фитнесс функции
      
#            self.RMS_S21_dB_DSA = RMS_S21_dB_DSA;  # Ррасчет СКО по амплитуде для ЦАТТ из листа List_Of_MF_S21_dB  -- Для фитнесс функции
#            self.RMS_S21_deg_DSA = RMS_S21_deg_DSA;  # Ррасчет СКО по градусам для ЦАТТ из листа List_Of_MF_S21_deg -- Для фитнесс функции

#            self.RMS_S21_dB_PhSh = RMS_S21_dB_PhSh;  # Расчет СКО по амплитуде для ФВ из листа List_Of_MF_S21_dB  -- Для фитнесс функции
#            self.RMS_S21_deg_PhSh = RMS_S21_deg_PhSh;   # Ррасчет СКО по градусам для ФВ из листа List_Of_MF_S21_deg -- Для фитнесс функции

#            self.FitnessValue_DSA = FitnessValue_DSA;  # Расчте фитнесс функции для ЦАТТ
#            self.FitnessValue_PhSh = FitnessValue_PhSh;  # Расчте фитнесс функции для ФВ

     


#def GetStateSystemDevice(Permutaion):

#    List_Of_States = [];  # лист непосредственно состояний

#    for i in range(number_of_states):  # по количеству состояний
#        TempState = GetState(i,Permutaion);
#        List_Of_States.append(TempState);

#    #####Наполнили лист 64 состояниями
#    List_Of_Max_S11_dB = []; # лист максимальных значений C11  для всех состояний
#    List_Of_Max_S22_dB = []; #лист максимальных значений C22  для всех состояний
#    List_Of_MF_P1dB = [];  # лист П1дБ значений для всех состояний
#    List_Of_MF_S21_dB = [];   # Лист с21 дБ для расчета СКО
#    List_Of_MF_S21_deg = [];  # Лист с21 град для расчета СКО


#    for i in range(number_of_states):
#        Current_Max_S11_dB = List_Of_States[i].S11_dB_Max;
#        Current_Max_S22_dB = List_Of_States[i].S22_dB_Max;
#        Current_P1dB = List_Of_States[i].P1dB;
#        Current_MF_S21_dB = List_Of_States[i].S21_dB_MF;
#        Current_MF_S21_deg = List_Of_States[i].S21_deg_MF;

#        List_Of_Max_S11_dB.append(Current_Max_S11_dB)
#        List_Of_Max_S22_dB.append(Current_Max_S22_dB);
#        List_Of_MF_P1dB.append(Current_P1dB);
#        List_Of_MF_S21_dB.append(Current_MF_S21_dB);
#        List_Of_MF_S21_deg.append(Current_MF_S21_deg);
#    ########################################################################      
#    ########################################################################
#    ########################################################################

#    Max_S11_dB = max(List_Of_Max_S11_dB); #Выбор макисмального значения из листа List_Of_Max_S11_dB  - Для фитнесс функции
#    Max_S22_dB = max(List_Of_Max_S22_dB);  #Выбор макисмального значения из листа List_Of_Max_S22_dB   - - Для фитнесс функции

#    #Min_P1dB = min(List_Of_MF_P1dB);  # В опорном Выбор минимального значения из листа List_Of_MF_P1dB;   -- Для фитнесс функции
#    Min_P1dB = List_Of_MF_P1dB[0];  # В опроном только

#    List_Of_IP1dB = List_Of_MF_P1dB;
#    List_of_OP1dB = [];
#    for i in range(number_of_states):
#        Current_OP1dB = List_Of_IP1dB[i]+((List_Of_MF_S21_dB[i])-1);
#        List_of_OP1dB.append(Current_OP1dB);
    
#    #######################
#    #######################
#    RMS_S21_dB_DSA = Get_RMS_S21_dB_DSA(List_Of_MF_S21_dB);  # Расчет СКО по амплитуде для ЦАТТ из листа List_Of_MF_S21_dB  -- Для фитнесс функции
#    RMS_S21_deg_DSA = Get_RMS_S21_deg_DSA(List_Of_MF_S21_deg); # Ррасчет СКО по градусам для ЦАТТ из листа List_Of_MF_S21_deg -- Для фитнесс функции

#    RMS_S21_dB_PhSh = Get_RMS_S21_dB_PhSh(List_Of_MF_S21_dB);
#    RMS_S21_deg_PhSh = Get_RMS_S21_deg_PhSh(List_Of_MF_S21_deg);

#    FitnessValue_DSA = Get_FitnessValue_DSA(Max_S11_dB,Max_S22_dB, RMS_S21_dB_DSA, RMS_S21_deg_DSA, Min_P1dB);
#    FitnessValue_PhSh = 0;
#    ########################
#    ########################
#    ########################
#    CombinationOrder = List_Of_States[0].CombinationOrder;
#    print('Permuation ready '+str(Permutaion));

#    return StateSystemDevice( CombinationOrder,   #Комбинация

#            List_Of_States,   # лист непосредственно состояний
#            List_Of_Max_S11_dB,  # лист максимальных значений для всех состояний
#            List_Of_Max_S22_dB, # лист максимальных значений для всех состояний
#            List_Of_MF_P1dB,   # лист П1дБ значений для всех состояний

#            List_Of_MF_S21_dB,   # Лист с21 дБ для расчета СКО
#            List_Of_MF_S21_deg,  # Лист с21 град для расчета СКО

#            Max_S11_dB,                  #Выбор макисмального значения из листа List_Of_Max_S11_dB  - Для фитнесс функции
#            Max_S22_dB,            #Выбор макисмального значения из листа List_Of_Max_S22_dB   - - Для фитнесс функции

#            Min_P1dB,               #Выбор минимального значения из листа List_Of_MF_P1dB;   -- Для фитнесс функции
      
#            RMS_S21_dB_DSA,  # Ррасчет СКО по амплитуде для ЦАТТ из листа List_Of_MF_S21_dB  -- Для фитнесс функции
#            RMS_S21_deg_DSA,  # Ррасчет СКО по градусам для ЦАТТ из листа List_Of_MF_S21_deg -- Для фитнесс функции

#            RMS_S21_dB_PhSh,  # Расчет СКО по амплитуде для ФВ из листа List_Of_MF_S21_dB  -- Для фитнесс функции
#            RMS_S21_deg_PhSh,   # Ррасчет СКО по градусам для ФВ из листа List_Of_MF_S21_deg -- Для фитнесс функции

#            FitnessValue_DSA,  # Расчте фитнесс функции для ЦАТТ
#            FitnessValue_PhSh)  # Расчте фитнесс функции для ФВ);





#def Get_RMS_S21_dB_DSA(List_Of_MF_S21_dB):
  
#    #считаем RMS амплитуды
#    Ref_List_S21 = [];  # нормированный лист ослабления
#    for state_is in range(number_of_states): 
#        Ref_List_S21.append(List_Of_MF_S21_dB[state_is]-List_Of_MF_S21_dB[0]);  # нормируем ослабление

#    ErrorList = [];
#    for i in range(number_of_states):     
#        ErrorList.append(Ref_List_S21[i] - NormasAttenuationList[i]);
#    #cреднее значение ошибки
#    MeanErrorS21dB = statistics.mean(ErrorList)    
#    #лист суммы для 
#    SumList=[];
#    for sum_i in range(number_of_states):
#        i_element = pow((ErrorList[sum_i]-MeanErrorS21dB),2);
#        SumList.append(i_element);
#    Sum_up = np.sum(SumList);
#    RMS_S21_dB_DSA = np.sqrt(Sum_up/(number_of_states-1));

#    return RMS_S21_dB_DSA;



#def Get_RMS_S21_deg_DSA(List_Of_MF_S21_deg):

#    MeanS21Deg = statistics.mean(List_Of_MF_S21_deg);
#    SumListDeg=[];
#    for i in range(number_of_states):
#        i_element = pow((List_Of_MF_S21_deg[i]-MeanS21Deg),2);
#        SumListDeg.append(i_element);
#    Sum_up_Deg = np.sum(SumListDeg);
#    RMS_S21_deg_DSA = np.sqrt(Sum_up_Deg/(number_of_states));
    
#    return RMS_S21_deg_DSA;

#def Get_RMS_S21_deg_PhSh(List_Of_MF_S21_deg):

#    #считаем RMS фазы
#    PhaseList = [];
#    for state_is in range(number_of_states):
#        PhaseList.append(List_Of_MF_S21_deg[state_is]-List_Of_MF_S21_deg[0]); # нормированный лист фазы
#    #Unwrap        
#    PhaseList = np.deg2rad(PhaseList);
#    PhaseList = np.unwrap(PhaseList);
#    PhaseList = np.rad2deg(PhaseList);
#    #Set Unwrap to StateList
#    #for unwrap_state in range(0,64,1):
#        #StateList[unwrap_state].StatePhase=PhaseList[unwrap_state];
#    #Делаем лист ошибки
#    ErrorList = [];
#    for i in range(0,64,1):     
#        ErrorList.append(PhaseList[i] - NormalPhaseList[i]);
#    #cреднее значение ошибки
#    MeanErrorPhase = statistics.mean(ErrorList)    
#    #лист суммы для 
#    SumList=[];
#    for sum_i in range(0,64,1):
#        i_element = pow((ErrorList[sum_i]-MeanErrorPhase),2);
#        SumList.append(i_element);
#    Sum_up = np.sum(SumList)
#    RMS_S21_deg_PhSh = np.sqrt(Sum_up/(len(SumList)-1));

#    return RMS_S21_deg_PhSh;



#def Get_RMS_S21_dB_PhSh(List_Of_MF_S21_dB):
  
#    S21List = List_Of_MF_S21_dB;

#    MeanS21 = statistics.mean(S21List);
#    SumList=[];
#    for i in range(0,64,1):
#        i_element = pow((S21List[i]-MeanS21),2);
#        SumList.append(i_element);
#    Sum_up = np.sum(SumList);
#    RMS_S21_dB_PhSh = np.sqrt(Sum_up/(len(SumList)-1));

#    return RMS_S21_dB_PhSh;
####################################################################################


#def Get_FitnessValue_DSA(Max_S11_dB,Max_S22_dB, RMS_S21_dB_DSA, RMS_S21_deg_DSA, Min_P1dB):
    
#    Max_S11_dB = Max_S11_dB;

#    if abs(S11_User)>abs(Max_S11_dB):
#        S11_Fit = (abs(S11_User)-abs(Max_S11_dB))*1.333333;
#    else:
#        S11_Fit = 0;

#    Max_S22_dB = Max_S22_dB;

#    if abs(S22_User)>abs(Max_S22_dB):
#        S22_Fit = (abs(S22_User)-abs(Max_S22_dB))*1.333333;
#    else:
#        S22_Fit = 0;

    
#    RMS_S21_deg_DSA = RMS_S21_deg_DSA;

#    if abs(RMS_S21_deg_DSA)>abs(RMS_S21_deg_DSA_User):
#        RMS_S21_deg_DSA_Fit = (abs(RMS_S21_deg_DSA)-abs(RMS_S21_deg_DSA_User))*6.66666666;
#    else:
#        RMS_S21_deg_DSA_Fit = 0;



#    RMS_S21_dB_DSA = RMS_S21_dB_DSA;

#    if abs(RMS_S21_dB_DSA)>abs(RMS_S21_dB_DSA_User):

#        RMS_S21_dB_DSA_Fit = (abs(RMS_S21_dB_DSA)-abs(RMS_S21_dB_DSA_User))*76.923;
#    else:
#        RMS_S21_dB_DSA_Fit = 0;


#    Min_P1dB_Fit = 0;

#    if abs(Min_P1dB_User)>abs(Min_P1dB):

#        Min_P1dB_Fit = (abs(Min_P1dB_User)-abs(Min_P1dB))*2;
#    else:
#        Min_P1dB_Fit = 0;
    

#    FitnessValue_DSA = S11_Fit #+ S22_Fit + RMS_S21_deg_DSA_Fit + RMS_S21_dB_DSA_Fit + Min_P1dB_Fit;

#    return FitnessValue_DSA;
     
##############################################################
##############################################################
##############################################################


##Combination = [3,1,4,5,6,7,2,8];  ## 8!Комбинация ферзей и перестановки
###Combination[0]  - C0
### h(C0) =  fitness - описывает качество перестановки - количество ферзей под ударом
##h_C0 = statistics.mean(Combination);
##T0 = 100;  # температура

### T0 = 100, C0,
###  k = 0;
###2. Ck+1 = через функцию перестановки,  def(Combination)
### [3 -- 6]  одна перестановка!!!
### Tk = T0+1 = a*Tk
### а = 0...1; 0.95
### dH = H(cK+1) - H(cK)
### if dh =0 break
### dH<0, new Combinationt = Ck+1  следущая перестановка
### dH>=0 расчитать веростноть p(dH) = exp(-dH/Tk) 
### c этой вероятностью принять за исходную комбинацию предыдущю и 
### c нее вернуть новую комбинацию. 1-ph() - взять новоую



#def GetStateAnnealing(StateNumber):  ### Меняет рандомно 2 позиции ячеек

#    CellList, CascadeList= GetStateData(StateNumber);  ## получили данные для состояния
#    Cascade = CascadeList;

#    Cascade_permutations = list(itertools.permutations(Cascade, len(BitList))); # 720 лист этих каскадов - 1 обьект листа это лист с 6 ячейками он же будет одинаково выдавать если состояние не менял??
#    Cell_Permutaitions = list(itertools.permutations(CellList, len(BitList)))  # 720 состояний ячеек, отсюда можно взять последовательность
    
#    Cascade_permutations[0]
#    Cascade_permutations[0]; # Первая перестановка, в ней поменять рандомно
    
#    CombinationOrder = [];
#    for i in range(len(CellList)):
#        CombinationOrder.append(Cell_Permutaitions[PermutaionNumber][i].CellValue);
#    Number = StateNumber; # номер состояния
    
#    Network = rf.cascade_list(Cascade_permutations[PermutaionNumber]);  # Берем все данные оттуда из собранного каска да

#    Frequency = Network.f.tolist(); # лист частот
#    S11_dB = [];   ### точки в по частотам
#    S22_dB = [];
#    S21_dB = [];
#    S21_deg = [];
#    for f in range(len(Frequency)):
#        S11_dB.append(float(Network.s11[str(Frequency[f])+'hz'].s_db[...]));
#        S21_dB.append(float(Network.s21[str(Frequency[f])+'hz'].s_db[...]));
#        S22_dB.append(float(Network.s22[str(Frequency[f])+'hz'].s_db[...]));
#        S21_deg.append(float(Network.s21[str(Frequency[f])+'hz'].s_deg[...]))
#        #print('waiting');
#    S11_dB_Max = max(S11_dB);
#    S22_dB_Max = max(S22_dB);
#    MF = math.ceil(len(Frequency)/2);  # центральная ччастотта
#    S21_dB_MF  = S21_dB[MF]; # выбор индекса на центральной частотте
#    S21_deg_MF = S21_deg[MF]
#    P1dB = GetP1dBCascade(CellList);  # обернул в функцию
#    return State(Number, CellList, CascadeList, 
#                  P1dB,  # точка сжатия для данного состояния при каскадном включении       
#                  S11_dB,  # лист значений
#                  S22_dB,  # /-/-/-/-/-/-/-
#                  S21_dB,  # /-/-/-/-/-/-/-
#                  S21_deg, # /-/-/-/-/-/-/-
#                  S11_dB_Max, # максимальное значение из списка частот
#                  S22_dB_Max, # максимальное значение из списка частот
#                  S21_dB_MF, # дБ на центральной частоте
#                  S21_deg_MF, # Градус на центральной частоте)
#                  CombinationOrder);




#############################################################
##################################################################

#def GetStateSystemDevice(Permutaion):

#    List_Of_States = [];  # лист непосредственно состояний

#    for i in range(number_of_states):  # по количеству состояний
#        TempState = GetState(i,Permutaion);
#        List_Of_States.append(TempState);

#    #####Наполнили лист 64 состояниями
#    List_Of_Max_S11_dB = []; # лист максимальных значений C11  для всех состояний
#    List_Of_Max_S22_dB = []; #лист максимальных значений C22  для всех состояний
#    List_Of_MF_P1dB = [];  # лист П1дБ значений для всех состояний
#    List_Of_MF_S21_dB = [];   # Лист с21 дБ для расчета СКО
#    List_Of_MF_S21_deg = [];  # Лист с21 град для расчета СКО


#    for i in range(number_of_states):
#        Current_Max_S11_dB = List_Of_States[i].S11_dB_Max;
#        Current_Max_S22_dB = List_Of_States[i].S22_dB_Max;
#        Current_P1dB = List_Of_States[i].P1dB;
#        Current_MF_S21_dB = List_Of_States[i].S21_dB_MF;
#        Current_MF_S21_deg = List_Of_States[i].S21_deg_MF;

#        List_Of_Max_S11_dB.append(Current_Max_S11_dB)
#        List_Of_Max_S22_dB.append(Current_Max_S22_dB);
#        List_Of_MF_P1dB.append(Current_P1dB);
#        List_Of_MF_S21_dB.append(Current_MF_S21_dB);
#        List_Of_MF_S21_deg.append(Current_MF_S21_deg);
#    ########################################################################      
#    ########################################################################
#    ########################################################################

#    Max_S11_dB = max(List_Of_Max_S11_dB); #Выбор макисмального значения из листа List_Of_Max_S11_dB  - Для фитнесс функции
#    Max_S22_dB = max(List_Of_Max_S22_dB);  #Выбор макисмального значения из листа List_Of_Max_S22_dB   - - Для фитнесс функции

#    #Min_P1dB = min(List_Of_MF_P1dB);  # В опорном Выбор минимального значения из листа List_Of_MF_P1dB;   -- Для фитнесс функции
#    Min_P1dB = List_Of_MF_P1dB[0];  # В опроном только

#    List_Of_IP1dB = List_Of_MF_P1dB;
#    List_of_OP1dB = [];
#    for i in range(number_of_states):
#        Current_OP1dB = List_Of_IP1dB[i]+((List_Of_MF_S21_dB[i])-1);
#        List_of_OP1dB.append(Current_OP1dB);
    
#    #######################
#    #######################
#    RMS_S21_dB_DSA = Get_RMS_S21_dB_DSA(List_Of_MF_S21_dB);  # Расчет СКО по амплитуде для ЦАТТ из листа List_Of_MF_S21_dB  -- Для фитнесс функции
#    RMS_S21_deg_DSA = Get_RMS_S21_deg_DSA(List_Of_MF_S21_deg); # Ррасчет СКО по градусам для ЦАТТ из листа List_Of_MF_S21_deg -- Для фитнесс функции

#    RMS_S21_dB_PhSh = Get_RMS_S21_dB_PhSh(List_Of_MF_S21_dB);
#    RMS_S21_deg_PhSh = Get_RMS_S21_deg_PhSh(List_Of_MF_S21_deg);

#    FitnessValue_DSA = Get_FitnessValue_DSA(Max_S11_dB,Max_S22_dB, RMS_S21_dB_DSA, RMS_S21_deg_DSA, Min_P1dB);
#    FitnessValue_PhSh = 0;
#    ########################
#    ########################
#    ########################
#    CombinationOrder = List_Of_States[0].CombinationOrder;
#    print('Permuation ready '+str(Permutaion));

#    return StateSystemDevice( CombinationOrder,   #Комбинация

#            List_Of_States,   # лист непосредственно состояний
#            List_Of_Max_S11_dB,  # лист максимальных значений для всех состояний
#            List_Of_Max_S22_dB, # лист максимальных значений для всех состояний
#            List_Of_MF_P1dB,   # лист П1дБ значений для всех состояний

#            List_Of_MF_S21_dB,   # Лист с21 дБ для расчета СКО
#            List_Of_MF_S21_deg,  # Лист с21 град для расчета СКО

#            Max_S11_dB,                  #Выбор макисмального значения из листа List_Of_Max_S11_dB  - Для фитнесс функции
#            Max_S22_dB,            #Выбор макисмального значения из листа List_Of_Max_S22_dB   - - Для фитнесс функции

#            Min_P1dB,               #Выбор минимального значения из листа List_Of_MF_P1dB;   -- Для фитнесс функции
      
#            RMS_S21_dB_DSA,  # Ррасчет СКО по амплитуде для ЦАТТ из листа List_Of_MF_S21_dB  -- Для фитнесс функции
#            RMS_S21_deg_DSA,  # Ррасчет СКО по градусам для ЦАТТ из листа List_Of_MF_S21_deg -- Для фитнесс функции

#            RMS_S21_dB_PhSh,  # Расчет СКО по амплитуде для ФВ из листа List_Of_MF_S21_dB  -- Для фитнесс функции
#            RMS_S21_deg_PhSh,   # Ррасчет СКО по градусам для ФВ из листа List_Of_MF_S21_deg -- Для фитнесс функции

#            FitnessValue_DSA,  # Расчте фитнесс функции для ЦАТТ
#            FitnessValue_PhSh)  # Расчте фитнесс функции для ФВ);

######################################################################



#BruteForce = [];
#for i in range(0,3,1):
#    BruteForce.append(GetStateSystemDevice(i))

#pickle.dump(BruteForce, open('DSA_BruteForce_1.pkl', 'wb'), protocol=pickle.HIGHEST_PROTOCOL);

#BruteForce =  pickle.load(open('BruteForce_S11.pkl', 'rb'));

#BruteForce = sorted(BruteForce, key = lambda StateSystemDevice: StateSystemDevice.Max_S11_dB);

#Sorted_Brute = [];

#Level = -7.2
#RMS = 0.54

#for i in range(len(BruteForce)):
#    if (BruteForce[i].Max_S11_dB<Level and BruteForce[i].Max_S22_dB<Level and BruteForce[i].RMS_S21_dB_DSA<RMS):
#        Sorted_Brute.append(BruteForce[i])
#    else:
#        continue;

#Sorted_Brute = sorted(Sorted_Brute, key = lambda StateSystemDevice: StateSystemDevice.RMS_S21_dB_DSA);

#for i in range(len(Sorted_Brute)):
#    print(Sorted_Brute[i].CombinationOrder)

#for i in range(len(Sorted_Brute)):
#    print(Sorted_Brute[i].RMS_S21_dB_DSA)
#print('sho')


#import xlsxwriter
#workbook = xlsxwriter.Workbook('Nik_Vuchich_Hello_For_DSA.xlsx');

#worksheet = workbook.add_worksheet();
#worksheet.write(0, 0, 'Сombination');
#worksheet.write(0, 1, 'Max_S11_dB');
#worksheet.write(0, 2, 'Max_S22_dB');
#worksheet.write(0, 3, 'RMS_S21_dB_DSA');
#worksheet.write(0, 4, 'RMS_S21_deg_DSA');
#worksheet.write(0, 5, 'RMS_S21_dB_PhSh');
#worksheet.write(0, 6, 'RMS_S21_deg_PhSh');

#for data in range(len(BruteForce)):  
#    row = data+1;   
#    worksheet.write(row, 0, str(BruteForce[data].CombinationOrder));
#    worksheet.write(row, 1, BruteForce[data].Max_S11_dB);
#    worksheet.write(row, 2, BruteForce[data].Max_S22_dB);
#    worksheet.write(row, 3, BruteForce[data].RMS_S21_dB_DSA);
#    worksheet.write(row, 4, BruteForce[data].RMS_S21_deg_DSA);
#    worksheet.write(row, 5, BruteForce[data].RMS_S21_dB_PhSh);
#    worksheet.write(row, 6, BruteForce[data].RMS_S21_deg_PhSh);
    
#workbook.close()
















    

#Combination = [3,1,4,5,6,7,2,8];  ## 8!Комбинация ферзей и перестановки
##Combination[0]  - C0
## h(C0) =  fitness - описывает качество перестановки - количество ферзей под ударом
#h_C0 = statistics.mean(Combination);

#T0 = 100;  # температура

## T0 = 100, C0,
##  k = 0;
##2. Ck+1 = через функцию перестановки,  def(Combination)
## [3 -- 6]  одна перестановка!!!
## Tk = T0+1 = a*Tk
## а = 0...1; 0.95
## dH = H(cK+1) - H(cK)
## if dh =0 break
## dH<0, new Combinationt = Ck+1  следущая перестановка
## dH>=0 расчитать веростноть p(dH) = exp(-dH/Tk) 
## c этой вероятностью принять за исходную комбинацию предыдущю и 
## c нее вернуть новую комбинацию. 1-ph() - взять новоую
