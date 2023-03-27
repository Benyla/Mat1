from sympy import *
from sympy.abc import x,y,z,v,u,w
from IPython.display import display, Math, Latex
init_printing()


def flux(V, r, firstint, secondint):
    if firstint[0] == symbols('u'):
        uint = firstint
        vint = secondint
    else:
        vint = firstint
        uint = secondint
    vruv = V.subs(x, r[0]).subs(y, r[1]).subs(z, r[2])
    nf = diff(r,u).cross(diff(r,v))
    return simplify(integrate(integrate(vruv.dot(nf), (u,uint[1],uint[2])), (v,vint[1],vint[2])))


def jacobi_kurve(r, var):
    rdiff = diff(r,var)
    if len(rdiff) == 2:
        return simplify(sqrt(rdiff[0]**2+rdiff[1]**2))
    else:
        return simplify(sqrt(rdiff[0]**2+rdiff[1]**2+rdiff[2]**2))


def flowkurve_w_bb(V, tid, rt):
    x = Function("x")
    y = Function("y")
    z = Function("z")
    t = symbols('t')
    eq1 = Eq(diff(x(t),t), V[0])
    eq2 = Eq(diff(y(t),t), V[1])
    eq3 = Eq(diff(z(t),t), V[2])
    sol = dsolve([eq1, eq2, eq3],ics={x(tid) : rt[0], y(tid) : rt[1], z(tid) : rt[2]})    
    return Matrix([sol[0].rhs, sol[1].rhs, sol[2].rhs])

def drej_om_z_akse(r):
    Q = Matrix([[cos(w), sin(w), 0], [-sin(w), cos(w), 0], [0,0,1]])
    return Q*r

def grad(f):
    return Matrix([diff(f,x), diff(f,y)])

def stat_punkter(grad):
    eq1 = Eq(grad[0], 0)
    eq2 = Eq(grad[1], 0)
    return solve([eq1, eq2], [x,y])

def partil_afledede(f, var, i):
    return diff(f, var, i)

def hesse_matrice(f):
    return Matrix([[diff(f, x, 2), diff(diff(f, x), y)], [diff(diff(f, y), x), diff(f, y, 2)]])


def analyse_af_stat_punkter(hesse_matrice, x0s, y0s):
    sol = []
    for i in range(len(x0s)):
        hesse_i_punkt = hesse_matrice.subs(x, x0s[i]).subs(y, y0s[i])
        sol.append(hesse_i_punkt.eigenvects()[0][0])
        sol.append(hesse_i_punkt.eigenvects()[1][0])
    for i in range(len(x0s)):
        if sol[i].evalf() > 0 and sol[i+1].evalf() > 0:
            print('Funktionsværdien i punktet: ('+str(x0s[i])+','+str(y0s[i])+') er et egentligt lokalt minimum')
        elif sol[i].evalf() < 0 and sol[i+1].evalf() < 0:
            print('Funktionsværdien i punktet: ('+str(x0s[i])+','+str(y0s[i])+') er et egentligt lokalt maksimum')
        elif (sol[i].evalf()*sol[i+1].evalf()) < 0:
            print('Funktionsværdien i punktet: ('+str(x0s[i])+','+str(y0s[i])+') er hverken lokal minimumsværdi eller maksimumsværdi')
        else:
            print('Ikke tilstrækkelig information. Kræver nærmere undersøgelse')



def grad_to_func(grad, funcval):
    K,c = symbols('K,c')
    function = (integrate(grad[0], x))+(integrate(grad[1]-(diff(integrate(grad[0], x),y)),y))
    if funcval == None:
        return function+K
    else:
        konstant = solve(Eq((function.subs(x, funcval[0]).subs(y, funcval[1]))+c, funcval[2]), c)
        final = function + konstant[0]
        return final


def approksimerende_andengrads_pol(f, udviklingspunkt):
    x0 = udviklingspunkt[0]
    y0 = udviklingspunkt[1]
    first = (f.subs(x, x0).subs(y, y0))+((diff(f,x).subs(x, x0).subs(y, y0))*(x-x0))+((diff(f,y).subs(x, x0).subs(y, y0))*(y-y0))
    second = (S(1)/2*(diff(f,x,2).subs(x, x0).subs(y, y0))*(x-x0)**2)
    third = (diff(diff(f,x),y).subs(x, x0).subs(y, y0))*(x-x0)*(y-y0)
    fourth = (S(1)/2*(diff(f,y,2).subs(x, x0).subs(y, y0))*(y-y0)**2)
    return first + second + third + fourth


def div(V):
    return diff(V[0], x) + diff(V[1], y) + diff(V[2], z)


def jacobi_rum(r):
    return simplify(det(Matrix([[diff(r[0], u), diff(r[0], v), diff(r[0], w)], 
                       [diff(r[1], u), diff(r[1], v), diff(r[1], w)], 
                       [diff(r[2], u), diff(r[2], v), diff(r[2], w)]])))

def kugle_rumfang(r):
    return (S(4)/3)*pi*r**3

def rumintegral_vol(r, uint, vint, wint):
    return integrate(jacobi_rum(r), (u,uint[1],uint[2]), (v,vint[1],vint[2]), (w,wint[1],wint[2]))


def fladeintegral_over_func(f, r, uint, vint, text=True):
    fruv = f.subs(x, r[0]).subs(y, r[1]).subs(z, r[2])
    sol = integrate(fruv*jacobi_flade(r), (u, uint[1], uint[2]), (v, vint[1], vint[2]))
    if text == None:
        return sol
    string1 = r'''Fladeintegralet \: af \: funktionen \: f\left(x,y,z\right) \: over \: den \: parametriserede \: falde \: F_{r} \: defineres \: ved:'''
    string2 = r'''\int_{F_{r}} f \: d\mu = \int_{c}^{d}\int_{a}^{b}f\left(r\left(u,v\right)\right)\mathit{Jacobi}_{r}\left(u,v\right)\mathit{dudv}'''
    string3 = r'''Her \: er \: \mathit{Jacobi}_{r}\left(u,v\right) \: givet \: ved:'''
    string4 = r'''\mathit{Jacobi}_{\boldsymbol{\mathit{r}}}\left(u,v\right)=| \: r_{u}^{\prime}\left(u,v\right)\times r_{v}^{\prime}\left(u,v\right)| '''
    string5 = r'''Jeg \: bestemmer \: først \: f\left(r\left(u,v\right)\right) \: og \: \mathit{Jacobi}_{\boldsymbol{\mathit{r}}}\left(u,v\right) '''
    string6 = r'''f\left(r\left(u,v\right)\right) = {%s}''' % (latex(fruv))
    string7 = r'''\mathit{Jacobi}_{\boldsymbol{\mathit{r}}}\left(u,v\right)=| \: r_{u}^{\prime}\left(u,v\right)\times r_{v}^{\prime}\left(u,v\right)| = {%s} ''' % (latex(jacobi_flade(r)))
    string8 = r'''Nu \: indsættes \: i \: formlen \: sammen \: med \: grænserne \: u \in [%s ,%s]\mathit{og} \: v \in [%s ,%s]''' % (latex(uint[1]), latex(uint[2]), latex(vint[1]), latex(vint[2]))
    string9 = r'''\int_{F_{r}} f \: d\mu = \int_{%s}^{%s}\int_{%s}^{%s}{%s}{%s}\mathit{dudv} = {%s}''' % (latex(uint[1]), latex(uint[2]), latex(vint[1]), latex(vint[2]), latex(fruv), latex(jacobi_flade(r)), latex(sol))
    display(Math(string1))
    display(Math(string2))
    display(Math(string3))
    display(Math(string4))
    display(Math(string5))
    display(Math(string6))
    display(Math(string7))
    display(Math(string8))
    display(Math(string9))
    return sol


def fladeintegral_areal(r, uint, vint):
    return integrate(jacobi_flade(r), (u, uint[1], uint[2]), (v, vint[1], vint[2]))

def jacobi_flade(r):
    return simplify((diff(r, u).cross(diff(r,v))).norm())

def tangentielle_kurveintegral(V, r, var, tstart, tslut):
    vru = V.subs(x, r[0]).subs(y, r[1]).subs(z, r[2])
    return integrate(vru*jacobi_kurve(r, var), (var, tstart, tslut))
