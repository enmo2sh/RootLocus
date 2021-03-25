import math
import sympy
from sympy import symbols, solve, diff, true, expand, simplify
import numpy
import matplotlib.pyplot as plt

def getRealPoles():
    realPoles=list()
    for p in range(0,len(poles)):
        if poles[p].imag==0 :
            realPoles.append(poles[p].real)
    realPoles.sort()
    for j in range(0,len(realPoles)-1):
        plt.plot([realPoles[j],realPoles[j+1]],[0,0],color='blue')
    return realPoles

def getalpha():
    sumOfPoles = 0
    sumOfZeros = 0
    for p in range(0, n):
        sumOfPoles = sumOfPoles + poles[p]
    for z in range(0, m):
        sumOfZeros = sumOfZeros + zeros[z]
    alpha = (sumOfPoles - sumOfZeros) / (n - m)
    print("a= ",alpha)
    return alpha

def getAsymptoteAngles():
    q = abs(n - m) - 1
    angles = list()
    for x in range(0, q + 1):
        temp = (180 * (2 * x + 1)) / (n - m)
        angles.append(temp);
    print("Asymptote Angles: ",angles)
    return angles


def breakAwayPoint():
    dk = diff(func)
    dk = simplify(dk)
    print("the derivative of K :", dk)
    root = solve(dk)
    for i in range(0,len(root)):
        root[i] = complex(root[i])
    print("roots :", root)
    realPoles=getRealPoles()
    result=list()
    for r in range(0,len(root)):
        for real in range (0,len(realPoles)-1):
             if root[r].real>realPoles[real]:
                if root[r].real<realPoles[real+1]:
                    result.append(root[r])
    print("Break Away Points :", result)
    return result


def getDepartureAngle(angle):
    Dangles = []
    count = 0
    for x in poles:
        angle = 0
        if x.imag != 0:
            for y in poles:
                if x != y:
                    if x.real == y.real:
                        angle = angle + 90
                    else:
                        temp = math.degrees(math.atan((x.imag - y.imag) / (x.real - y.real)))
                        if x.real > y.real:
                            if x.imag>y.imag:
                                angle = angle + temp
                            else:
                                angle=angle+360+temp
                        else:
                            if x.imag > y.imag:
                                angle = angle + 180 + temp
                            else:
                                angle = angle + 180 + temp
                plt.plot([x.real, y.real], [x.imag, y.imag],color='green')
            angle=180-angle
            if angle>360:
                angle=angle-360
            elif angle<0:
                angle=angle+360
            if(x.imag<0):
                angle=angle+180
            Dangles.append([])
            Dangles[count].append(x)
            Dangles[count].append(angle)
            count = count + 1
    print("Departure Angles " , Dangles)
    return (Dangles)

def getIntersectionWithImagAxis():
    coeffs=sympy.Poly(func).all_coeffs()
    k=sympy.Symbol('k')
    coeffs[len(coeffs)-1]=coeffs[len(coeffs)-1]+k
    routhTable=[]
    routhTable.append([])
    routhTable.append([])
    count=0
    while(count<len(coeffs)-1):
        routhTable[0].append(coeffs[count])
        routhTable.append([])
        routhTable[1].append(coeffs[count+1])
        count=count+2
        if(count==len(coeffs)-1):
            routhTable[0].append(coeffs[count])
    routhTable[0].append(0)
    while(len(routhTable[1])<len(routhTable[0])):
        routhTable[1].append(0)
    for r in range (2,len(coeffs)):
        routhTable.append([])
        for e in range(0,len(routhTable[r-1])-1):
            term1=(routhTable[r-1][e]*routhTable[r-2][e+1])
            term2=(routhTable[r-1][e+1]*routhTable[r-2][e])
            if(term1-term2)==0:
                routhTable[r].append(0)
            else:
                routhTable[r].append((term1-term2)/routhTable[r-1][e])
    print("Routh Table" , routhTable)
    equationOfk = routhTable[len(coeffs)-2][0]
    value = sympy.solve(equationOfk)
    equ = routhTable[len(coeffs) - 3][0] * s ** 2
    equ =equ+ routhTable[len(coeffs) - 3][1]
    equ2=equ.subs(k,value[0])
    intersection=solve(equ2)
    for i in range(0,len(intersection)):
        intersection[i]=complex(intersection[i])
    print("intersection with imag axis :" , intersection)
    return intersection

def Draw():
    x = list()
    y = list()
    for p in range(0, n):
        x.append(poles[p].real)
        y.append(poles[p].imag)
    plt.scatter(x, y, label="stars", color="green",
                marker="*", s=50)
    for a in range(0, len(angles)):
        endy = 0 + 125 * numpy.math.sin(math.radians(angles[a]))
        endx = alpha + 125 * math.cos(math.radians(angles[a]))
        plt.plot([alpha, endx], [0, endy],color='red',ls = '--')
    for a in range(0,len(Dangles)):
        yf = intersection[1-a].imag
        plt.plot(intersection[1-a].real, intersection[1-a].imag+(1.5*(-1*a)), 'o', color="red")
        if(yf>0):
            xf=Dangles[a][0].real-((yf-Dangles[a][0].imag)/math.tan(math.radians(180-Dangles[a][1])))
        else:
            xf = Dangles[a][0].real+((yf - Dangles[a][0].imag) / math.tan(math.radians(Dangles[a][1]-180)))
        plt.plot([Dangles[a][0].real,xf],[Dangles[a][0].imag,yf],color='blue')
        plt.plot([xf, xf + 50 * math.cos(math.radians(angles[a+1]))],
                 [yf, yf + 50 * numpy.math.sin(math.radians(angles[a+1]))],color='blue')
    plt.plot(breakPoints[0].real,0, 'o', color="orange")
    yp=numpy.linspace(-75,75,30)
    temp = 1 + (yp ** 2 / (alpha - breakPoints[0].real) ** 2)
    temp = temp * (alpha - breakPoints[0].real) ** 2
    xp = numpy.sqrt(temp) + alpha
    plt.plot(xp,yp,color='blue')
    plt.xlabel('real- axis')
    plt.ylabel('imag. - axis')
    plt.title('ROOT LOCUS')
    plt.grid(True)
    plt.show()



n = int(input("enter the number of poles: "))
m = int(input("enter the number of zeros: "))
poles = list()
zeros = list()
if (n > 0):
    print("enter the poles")
    for p in range(0, n):
        poles.append(complex(input()))
if (m > 0):
    print("enter the zeros")
    for z in range(0, m):
        zeros.append(complex(input()))
s = symbols('s')
func = (s - poles[0])
for g in range(1, n):
    func = (s - poles[g]) * func
func = expand(func)
print("Equation :" , func)
alpha = getalpha()
angles = getAsymptoteAngles()
breakPoints=breakAwayPoint()
Dangles=getDepartureAngle(0)
intersection=getIntersectionWithImagAxis()
Draw()

