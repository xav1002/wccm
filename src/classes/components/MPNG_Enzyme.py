import json

class MPNG_Enzyme:


    def __init__(self,entry,names,sysname,reactions,substrates,products):
        self.__entry = entry
        self.__names = names
        self.__sysname = sysname
        self.__reactions = reactions
        self.__substrates = substrates
        self.__products = products

        self.__rxn_reversibility = {}

    @property
    def entry(self) -> str:
        return self.__entry

    @property
    def names(self) -> list[str]:
        return self.__names

    @property
    def sysname(self) -> str:
        return self.__sysname

    @property
    def reactions(self) -> list[str]:
        return self.__reactions

    @property
    def substrates(self) -> list[str]:
        return self.__substrates

    @property
    def products(self) -> list[str]:
        return self.__products

    def set_rxn_reversibility(self,rxn_entry:str,rev:bool) -> None:
        self.__rxn_reversibility[rxn_entry] = rev

    def toJSON(self):
        return json.dumps({
                    'entry':self.__entry,
                    'names':self.__names,
                    'sysname':self.__sysname,
                    'reactions':self.__reactions,
                    'substrates':self.__substrates,
                    'products':self.__products,
                    # 'conc':self.__conc,
                    # 'explored':self.__explored
                })

    def fromJSON(dict_from_json):
        rxn = MPNG_Enzyme(dict_from_json['entry'],dict_from_json['names'],
                               dict_from_json['sysname'],dict_from_json['reactions'],
                               dict_from_json['substrates'],dict_from_json['products'])
        return rxn