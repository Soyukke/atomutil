class MaterialColor:
    def __init__(self, ambient, diffuse, specular, shininess):
        self.ambient = ambient
        self.diffuse = diffuse
        self.specular = specular
        self.shininess = shininess

def name2color(name):
    if name=="C":
        return obsidian
    elif name=="Si":
        return black_plastic       
    elif name=="H":
        return white_plastic
    elif name=="O":
        return ruby
    elif name=="N":
        return cyan_plastic
    elif name=="P":
        return cyan_rubber
    elif name=="S":
        return yellow_plastic
    elif name=="F":
        return green_plastic
    elif name=="Cl":
        return green_rubber
    elif name=="Br":
        return chrome
    elif name=="I":
        return pearl
    elif name=="Cu":
        return copper
    elif name=="Ag":
        return silver
    elif name=="Au":
        return gold
    elif name=="B":
        return brass
    elif name=="Al":
        return bronze
    else:
        return silver


#  ruby(ルビー)
ruby = MaterialColor([0.1745, 0.01175, 0.01175, 1.0],[0.61424, 0.04136, 0.04136, 1.0],[0.727811, 0.626959, 0.626959, 1.0],76.8)
#  emerald(エメラルド)
emerald = MaterialColor([0.0215, 0.1745, 0.0215, 1.0],
[0.07568, 0.61424, 0.07568, 1.0],
[0.633, 0.727811, 0.633, 1.0],
    76.8)
    #  jade(翡翠)
jade = MaterialColor([0.135, 0.2225, 0.1575, 1.0],
[0.54, 0.89, 0.63, 1.0],
[0.316228, 0.316228, 0.316228, 1.0],
    12.8)
    #  obsidian(黒曜石)
obsidian = MaterialColor([0.05375, 0.05, 0.06625, 1.0],
[0.18275, 0.17, 0.22525, 1.0],
[0.332741, 0.328634, 0.346435, 1.0],
    38.4)
    #  pearl(真珠)
pearl = MaterialColor([0.25, 0.20725, 0.20725, 1.0],
[1, 0.829, 0.829, 1.0],
[0.296648, 0.296648, 0.296648, 1.0],
    10.24)#  turquoise(トルコ石)
turquoise = MaterialColor([0.1, 0.18725, 0.1745, 1.0],
[0.396, 0.74151, 0.69102, 1.0],
[0.297254, 0.30829, 0.306678, 1.0],
    12.8)#  brass(真鍮)
brass = MaterialColor([0.329412, 0.223529, 0.027451, 1.0],
[0.780392, 0.568627, 0.113725, 1.0],
[0.992157, 0.941176, 0.807843, 1.0],
    27.89743616)#  bronze(青銅)
bronze = MaterialColor([0.2125, 0.1275, 0.054, 1.0],
[0.714, 0.4284, 0.18144, 1.0],
[0.393548, 0.271906, 0.166721, 1.0],
    25.6)#  chrome(クローム)
chrome = MaterialColor([0.25, 0.25, 0.25, 1.0],
[0.4, 0.4, 0.4, 1.0],
[0.774597, 0.774597, 0.774597, 1.0],
    76.8)#  copper(銅)
copper = MaterialColor([0.19125, 0.0735, 0.0225, 1.0],
[0.7038, 0.27048, 0.0828, 1.0],
[0.256777, 0.137622, 0.086014, 1.0],
    12.8)#  gold(金)
gold = MaterialColor([0.24725, 0.1995, 0.0745, 1.0],
[0.75164, 0.60648, 0.22648, 1.0],
[0.628281, 0.555802, 0.366065, 1.0],
    51.2)#  silver(銀)
silver = MaterialColor([0.19225, 0.19225, 0.19225, 1.0],
[0.50754, 0.50754, 0.50754, 1.0],
[0.508273, 0.508273, 0.508273, 1.0],
    51.2)#  プラスチック(黒)
black_plastic = MaterialColor([0.0, 0.0, 0.0, 1.0],
[0.01, 0.01, 0.01, 1.0],
[0.50, 0.50, 0.50, 1.0],
    32)#  プラスチック(シアン)
cyan_plastic = MaterialColor([0.0, 0.1, 0.06, 1.0],
[0.0, 0.50980392, 0.50980392, 1.0],
[0.50196078, 0.50196078, 0.50196078, 1.0],
    32)#  プラスチック(緑)
green_plastic = MaterialColor([0.0, 0.0, 0.0, 1.0],
[0.1, 0.35, 0.1, 1.0],
[0.45, 0.55, 0.45, 1.0],
    32)#  プラスチック(赤)
red_plastic = MaterialColor([0.0, 0.0, 0.0, 1.0],
[0.5, 0.0, 0.0, 1.0],
[0.7, 0.6, 0.6, 1.0],
    32)#  プラスチック(白)
white_plastic = MaterialColor([0.0, 0.0, 0.0, 1.0],
[0.55, 0.55, 0.55, 1.0],
[0.70, 0.70, 0.70, 1.0],
    32)#  プラスチック(黄)
yellow_plastic = MaterialColor([0.0, 0.0, 0.0, 1.0],
[0.5, 0.5, 0.0, 1.0],
[0.60, 0.60, 0.50, 1.0],
    32)#  ゴム(黒)
black_rubber = MaterialColor([0.02, 0.02, 0.02, 1.0],
[0.01, 0.01, 0.01, 1.0],
[0.4, 0.4, 0.4, 1.0],
    10.0)#  ゴム(シアン)
cyan_rubber = MaterialColor([0.0, 0.05, 0.05, 1.0],
[0.4, 0.5, 0.5, 1.0],
[0.04, 0.7, 0.7, 1.0],
    10.0)#  ゴム(緑)
green_rubber = MaterialColor([0.0, 0.05, 0.0, 1.0],
[0.4, 0.5, 0.4, 1.0],
[0.04, 0.7, 0.04, 1.0],
    10.0)#  ゴム(赤)
red_rubber = MaterialColor([0.05, 0.0, 0.0, 1.0],
[0.5, 0.4, 0.4, 1.0],
[0.7, 0.04, 0.04, 1.0],
    10.0)#  ゴム(白)
white_rubber = MaterialColor([0.05, 0.05, 0.05, 1.0],
[0.5, 0.5, 0.5, 1.0],
[0.7, 0.7, 0.7, 1.0],
    10.0)#  ゴム(黄)
yellow_rubber = MaterialColor([0.05, 0.05, 0.0, 1.0],
[0.5, 0.5, 0.4, 1.0],
[0.7, 0.7, 0.04, 1.0],
    10.0)
