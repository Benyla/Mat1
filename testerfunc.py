from sympy import *
from sympy.abc import v,u,w,x,y,z
from IPython.display import display, Math, Latex
import inspect
init_printing()




# non math related

#Hvis den bruges i en funktion, der kaldes inden i en anden funktion fucker det op
def retrieve_name(var):
    callers_local_vars = inspect.currentframe().f_back.f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]

def num_of_vars(f):
    vars = list(f.free_symbols)
    return len(vars)




#Fluxen gennem flade
def flux_gennem_flade(V, r, firstint, secondint):
    if firstint[0] == symbols('u'):
        uint = firstint
        vint = secondint
    else:
        vint = firstint
        uint = secondint
    vruv = V.subs(x, r[0]).subs(y, r[1]).subs(z, r[2])
    nf = diff(r,u).cross(diff(r,v))
    return simplify(integrate(integrate(vruv.dot(nf), (u,uint[1],uint[2])), (v,vint[1],vint[2])))

#Jacobi-funktion for en kruve
def jacobi_kurve(r, u, text=True):
        sol = simplify((Matrix([diff(r[0],u), diff(r[1],u), diff(r[2],u)])).norm())
        if text==False:
                return sol
        string1 = r'''Jacobi-funktionen \: for \: en \: kurve \: i \: planen \: og \: rummet \: er \: givet \: ved: \mathit{Jacobi}_{r}\left(u\right)=|r'\left(u\right)| = |\left[\begin{array}{c}
{%s}  
\\
 {%s}   
\\
 {%s}   
\end{array}\right]| = {%s}''' % (latex(diff(r[0],u)), latex(diff(r[1],u)), latex(diff(r[2],u)), latex(sol))
        display(Math(string1))
        return sol

#fiks, så kan bruges uden med x frem for x(t)
def flowkurve_w_bb_rum(V, tid, rt):
    t = symbols('t')
    eq1 = Eq(diff(x(t),t), V[0])
    eq2 = Eq(diff(y(t),t), V[1])
    eq3 = Eq(diff(z(t),t), V[2])
    sol = dsolve([eq1, eq2, eq3],ics={x(tid) : rt[0], y(tid) : rt[1], z(tid) : rt[2]})    
    return Matrix([sol[0].rhs, sol[1].rhs, sol[2].rhs])

def flowkurve_w_bb_plan(V, tid, rt):
    t = symbols('t')
    eq1 = Eq(diff(x(t),t), V[0])
    eq2 = Eq(diff(y(t),t), V[1])
    sol = dsolve([eq1, eq2],ics={x(tid) : rt[0], y(tid) : rt[1]})    
    return Matrix([sol[0].rhs, sol[1].rhs])

def drej_om_z_akse(r):
    Q = Matrix([[cos(w), sin(w), 0], [-sin(w), cos(w), 0], [0,0,1]])
    return Q*r


#Tager en funktion som input. Kan kun bruges til funktioner af x,y og z for nu.
def grad(f):
    number_of_vars = num_of_vars(f)
    var_list = [x, y, z]
    sol = Matrix([])
    for i in range(number_of_vars):
        sol = sol.row_insert(sol.shape[0], Matrix([diff(f, var_list[i])]))
    return sol

#Simpel niveaukurve. Tager en funktion og en konstant, som indput. 
def niveaukurve(f, c, text=True):
    if text == False:
        return simplify(Eq(f,c))
    string1 = r'''Niveaukurven \: for \: {%s} \: svarende \: til \: niveauet \: {%s} \: er \: givet \: ved: {%s}''' % (retrieve_name(f)[0], latex(c), latex(simplify(Eq(f,c))))
    display(Math(string1))

    return simplify(Eq(f,c))

#Funktionen her forsøger at give et udtryk for niveaukurven, som funktion af både x og y. 
#Tager en funktion og konstant som input.
def niveaukurve_undersøgelse(f, c):
    try:
        eq1 = simplify(Eq(solve(niveaukurve(f, c, text=False),y)[0], y))
        string1 = r'''Niveaukurven \: for \: {%s} \: svarende \: til \: niveauet \: {%s} \: som \: funktion \: af \: x \: er \: givet \: ved: {%s}''' % (retrieve_name(f)[0], latex(c), latex(eq1))
        display(Math(string1))
    except:
        print('NaN')
    try:
        eq2 = simplify(Eq(solve(niveaukurve(f, c, text=False),x)[0], x))
        string2 = r'''Niveaukurven \: for \: {%s} \: svarende \: til \: niveauet \: {%s} \: som \: funktion \: af \: y \: er \: givet \: ved: {%s}''' % (retrieve_name(f)[0], latex(c), latex(eq2))
        display(Math(string2))
    except:
        print('Nan')


#grad, retningsvektor og punkt matrix eller tuple. funcname skal defineres, ex: func = 'f'
def retningsaflede_i_punkt(grad, retningsvektor, punkt, funcname, text=True):    
    var_list = [x,y,z]
    for i in range(len(punkt)):
        if i < 1:
            grad_i_punkt = grad.subs(var_list[i], punkt[i])
        else:
            grad_i_punkt = grad_i_punkt.subs(var_list[i], punkt[i])
    enhedsretningsvektor = retningsvektor/retningsvektor.norm()

    if text == False:
        return simplify(grad_i_punkt.dot(enhedsretningsvektor))
    
    string1 = r'''Den \: retningsaflede \: af \: {%s} \: i \: punktet \: {%s} \: i \: retningen \: givet \: ved \: retningsvektoren \: v = {%s} \: er \: givet \: ved: {%s}''' % (funcname, latex(punkt), latex(retningsvektor), latex(simplify(grad_i_punkt.dot(enhedsretningsvektor))))
    display(Math(string1))

    return simplify(grad_i_punkt.dot(enhedsretningsvektor))

#grad og punkt matrix eller tuple. funcname skal defineres, ex: func = 'f'
#Finder maksimale og minimale værdi af den retningsaflede i et givent punkt
def bestem_maks_min_af_retningsaflede_i_punkt(grad, punkt, funcname, text=True):
    var_list = [x,y,z]
    for i in range(len(punkt)):
        if i < 1:
            grad_i_punkt = grad.subs(var_list[i], punkt[i])
        else:
            grad_i_punkt = grad_i_punkt.subs(var_list[i], punkt[i])
    enhedsretningsvektor = Matrix([cos(u), sin(u)])
    expr = grad_i_punkt.dot(enhedsretningsvektor)
    sol = [s for s in solveset(diff(expr, u), u, Interval(-pi, pi))]
    maks_min_liste = []
    for i in range(len(sol)):
        retning = Matrix([cos(sol[i]), sin(sol[i])])
        maks_min_liste.append(grad_i_punkt.dot(retning))   

    if text == False:
        return maks_min_liste
    
    string1 = r'''Lad \: e \: være \: en \: retningsvektor \: givet \: ved \: e =(\cos (u),\sin (u))'''
    string2 = r'''Gradienten \: i \: punktet \: er \: givet \: ved \: \nabla {%s} \left(x_{o},y_{0}\right) = \nabla {%s} \left({%s},{%s}\right) = {%s}''' % (funcname, funcname, latex(punkt[0]), latex(punkt[1]), latex(grad_i_punkt))
    string3 = r'''Vi \: kan \: da \: finde \: et \: udtryk \: for \: den \: retningsaflede \: i \: alle \: retninger \: ved: \nabla {%s} \left({%s},{%s}\right) \cdot e = {%s}''' %(funcname,latex(punkt[0]), latex(punkt[1]), latex(expr))
    string4 = r'''Jeg \: finder \: nu \: løsningerne \: for \: \frac{\partial}{\partial u}\left({%s}\right) = 0 \leftrightarrow u = {%s} \: i \: intervallet \: [-\pi,\pi] ''' % (latex(expr), latex(sol))
    string5 = r'''Nu \: kan \: løsningerne \: indsættes \: i \: e \: og \: derefter \: findes \: prikproduktet \: mellem \: \nabla {%s} \left({%s},{%s}\right) \: og \: e''' % (funcname, latex(punkt[0]), latex(punkt[1]))
    string6 = r'''Den \: maksiamel \: værdi \: af \: den \: retningsaflede \: i \: punktet \: \left({%s},{%s}\right) \: er \: altså \: givet \: ved \: {%s}''' %(latex(punkt[0]), latex(punkt[1]), latex(max(maks_min_liste)))
    string7 = r'''Den \: mindste \: værdi \: af \: den \: retningsaflede \: i \: punktet \: \left({%s},{%s}\right) \: er \: givet \: ved \: {%s}''' %(latex(punkt[0]), latex(punkt[1]), latex(min(maks_min_liste)))
    display(Math(string1))
    display(Math(string2))
    display(Math(string3))
    display(Math(string4))
    display(Math(string5))
    display(Math(string6))
    display(Math(string7))

    return maks_min_liste


#Tager en gradient på matrix form. funcname er funktionens navn ex: funcname = 'f'
def stat_punkter_2_var(grad, funcname, text=True):
    eq1 = Eq(grad[0], 0)
    eq2 = Eq(grad[1], 0)
    if text == False:
        return solve([eq1, eq2], [x,y])
    string1 = r'''De \: stationære \: punkter \: findes \: ved at \: løse \: følgende \: ligningssystem:'''
    string2 = r'''\frac{\partial}{\partial x}{%s}=0''' % (funcname)
    string3 = r'''\frac{\partial}{\partial y}{%s}=0''' % (funcname)
    string4 = r'''Herved \: opstilles \: ligningssystemet:'''
    string5 = r'''{%s}=0''' % (latex(grad[0]))
    string6 = r'''{%s}=0''' % (latex(grad[1]))
    string7 = r'''Løsningerne \: udregnes \: til: {%s}''' % (latex(solve([eq1, eq2], [x,y])))
    string8 = r'''Følgende \: er \: altså \: stationære \: punkter \: for \: funktionen, {%s}:''' % (funcname)
    display(Math(string1))
    display(Math(string2))
    display(Math(string3))
    display(Math(string4))
    display(Math(string5))
    display(Math(string6))
    display(Math(string7))
    display(Math(string8))
    for i in range(len(solve([eq1, eq2], [x,y]))):
        stringx = r'''Punkt \: {%s} = {%s}''' % (latex(i), latex((solve([eq1, eq2], [x,y])[i])))
        display(Math(stringx))
    return solve([eq1, eq2], [x,y])


#Tager en gradient på matrix form. funcname er funktionens navn ex: funcname = 'f'
def stat_punkter_3_var(grad, funcname, text=True):
    eq1 = Eq(grad[0], 0)
    eq2 = Eq(grad[1], 0)
    eq3 = Eq(grad[2], 0)
    if text == False:
        return solve([eq1, eq2, eq3], [x,y,z])
    string1 = r'''De \: stationære \: punkter \: findes \: ved at \: løse \: følgende \: ligningssystem:'''
    string2 = r'''\frac{\partial}{\partial x}{%s}=0''' % (funcname)
    string3 = r'''\frac{\partial}{\partial y}{%s}=0''' % (funcname)
    string4 = r'''\frac{\partial}{\partial z}{%s}=0''' % (funcname)
    string5 = r'''Herved \: opstilles \: ligningssystemet:'''
    string6 = r'''{%s}=0''' % (latex(grad[0]))
    string7 = r'''{%s}=0''' % (latex(grad[1]))
    string8 = r'''{%s}=0''' % (latex(grad[2]))
    string9 = r'''Løsningerne \: udregnes \: til: {%s}''' % (latex(solve([eq1, eq2, eq3], [x,y,z])))
    string10 = r'''Følgende \: er \: altså \: stationære \: punkter \: for \: funktionen, {%s}:''' % (funcname)
    display(Math(string1))
    display(Math(string2))
    display(Math(string3))
    display(Math(string4))
    display(Math(string5))
    display(Math(string6))
    display(Math(string7))
    display(Math(string8))
    display(Math(string9))
    display(Math(string10))
    for i in range(len(solve([eq1, eq2, eq3], [x,y,z]))):
        stringx = r'''Punkt \: {%s} = {%s}''' % (latex(i), latex((solve([eq1, eq2, eq3], [x,y,z])[i])))
        display(Math(stringx))
    return solve([eq1, eq2, eq3], [x,y,z])

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


def divv(V, vars):
    dV = Matrix([diff(V[i], vars[i]) for i in range(len(V))])
    return sum(dV[i] for i in range(len(V)))



def jacobi_rum(r, u, v, w, text=True):
    sol = simplify(abs(det(Matrix([[diff(r[0], u), diff(r[0], v), diff(r[0], w)], 
                       [diff(r[1], u), diff(r[1], v), diff(r[1], w)], 
                       [diff(r[2], u), diff(r[2], v), diff(r[2], w)]]))))
    if text == False:
        return sol
    string1 = r'''Givet \: en \: parameterfremstilling, \: r, \: for \: et \: rumligt \: område, \: \Omega, \: i \: rummet. \: Den \: til \: r \: hørende \: Jacobifunktion \: udregnes \: således:'''
    string2 = r'''Jacobi_{r} = |\det\left(\left[\begin{array}{cc}
\frac{\partial x}{\partial u}r\left(u,v,w\right) & \frac{\partial x}{\partial v}r\left(u,v,w\right) & \frac{\partial x}{\partial w}r\left(u,v,w\right) 
\\
 \frac{\partial y}{\partial u}r\left(u,v,w\right) & \frac{\partial y}{\partial v}r\left(u,v,w\right) & \frac{\partial y}{\partial w}r\left(u,v,w\right)
\\
 \frac{\partial z}{\partial u}r\left(u,v,w\right) & \frac{\partial z}{\partial v}r\left(u,v,w\right) & \frac{\partial z}{\partial w}r\left(u,v,w\right)
\end{array}\right]\right)| = |\left[\begin{array}{ccc}
{%s}  & {%s}  & {%s}  
\\
 {%s}  & {%s}  & {%s}  
\\
 {%s}  & {%s}  & {%s}  
\end{array}\right]| = {%s}''' % (latex(diff(r[0], u)), latex(diff(r[0], v)), latex(diff(r[0], w)), latex(diff(r[1], u)), latex(diff(r[1], v)), latex(diff(r[1], w)), latex(diff(r[2], u)), latex(diff(r[2], v)), latex(diff(r[2], w)), latex(sol))
    display(Math(string1))
    display(Math(string2))
    return sol

def kugle_rumfang(r):
    return (S(4)/3)*pi*r**3

def rumintegral_vol(r, uint, vint, wint, text=True):
    sol = integrate(jacobi_rum(r, u, v, w, text=False), (u,uint[1],uint[2]), (v,vint[1],vint[2]), (w,wint[1],wint[2]))
    if text==False:
        return sol
    string1 = r'''Volumenet \: af \: den \: parametriserede \: rumlige \: område \: \Omega_{r} : r\left(u,v,w\right)=\left(x\left(u,v,w\right),y\left(u,v,w\right),z\left(u,v,w\right)\right)u\in \left[a,b\right] v\in \left[c,d\right] w\in \left[h,l\right]'''
    string2 = r'''defineres \: som \: rumintegralet: \int_{\Omega_{r}}^{}1d \mu = \int_{{h}}^{l}\int_{{c}}^{d}\int_{{a}}^{b}Jacobi_{r} \:d u d v dw'''
    display(Math(string1))
    display(Math(string2))
    string3 = r'''Hermed \: er \: Jacobi \: funktionen \: givet \: ved: {%s}''' % (latex(jacobi_rum(r, u, v, w)))
    display(Math(string3))
    string4 = r'''\int_{\Omega_{r}}^{}1d \mu = \int_{{{%s}}}^{{%s}}\int_{{{%s}}}^{{%s}}\int_{{{%s}}}^{{%s}}{%s}d u d v dw = {%s}''' % (latex(wint[1]), latex(wint[2]), latex(vint[1]), latex(vint[2]), latex(uint[1]), latex(uint[2]), latex(jacobi_rum(r, u, v, w, text=False)), latex(sol))
    display(Math(string4))
    return sol

#Bruges til at beregne masse, hvor funktionen f, er en masse tæthedsfunktion
def fladeintegral_af_func_over_flade(f, r, uint, vint, text=True):
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

#Tager en parameterfremstilling og 2 symboler.
def jacobi_flade(r,u,v, text=True):
    sol = simplify((diff(r, u).cross(diff(r, v))).norm())
    if text == False:
        return sol
    string1 = r'''Givet \: en \: parameterfremstilling, \: r, \: for \: en \: flade, \: F, \: i \: rummet. \: Den \: til \: r \: hørende \: Jacobifunktion \: udregnes \: således:'''
    string2 = r'''Jacobi_{r} = |\frac{\partial}{\partial u}r\left(u,v\right)\times \frac{\partial}{\partial v}r\left(u,v\right)| = |{%s} \times {%s}| = {%s}''' % (latex(diff(r, u)), latex(diff(r, v)), latex(sol))
    display(Math(string1))
    display(Math(string2))
    return sol

#var er var i parameterfremstilligen
def kurveintegral_af_f_over_kurve(f, r, var, start, slut):
    vru = f.subs(x, r[0]).subs(y, r[1]).subs(z, r[2])
    return integrate(vru*jacobi_kurve(r, var), (var, start, slut))

def rott(V):
    return Matrix([(diff(V[2], y)-diff(V[1], z)), (diff(V[0], z)-diff(V[2], x)), (diff(V[1], x)-diff(V[0], y))])

def tjek_v_gradfelt(V):
    if rott(V).norm() == 0:
        return True
    else:
        return False
    
def stam_til_vektorfelt(V):
    Vuxuyuz = Matrix([u*x,u*y,u*z])
    delres = Matrix([Vuxuyuz[0].subs(x, V[0]), Vuxuyuz[1].subs(y, V[1]), Vuxuyuz[2].subs(z, V[2])])
    sol = Matrix([x,y,z]).dot(Matrix([integrate(delres[0], (u,0,1)), integrate(delres[1], (u,0,1)), integrate(delres[2], (u,0,1))]))
    return sol


#Tager en matrix som input
def diagonalisering(A, text=True):
    data = A.eigenvects()
    geo_sum = 0
    for i in range(len(A.row(0))):
        geo_sum += data[i][1]
    if len(data) == len(A.row(0)) or geo_sum == len(A.row(0)):
        if text == True:
            print('kriteriet er opfyldt, for at der kan diagonaliseres ved similartransformation')
    else:
        if text == True:
            print('kriteriet for at der kan diagonaliseres ved similartransformation, er ikke opfyldt')
    
    egenvals = []
    Q_matrix = Matrix([])
    for i in range(len(A.row(0))):
        egenvals.append(data[i][0])
        divider = ((data[i][2][0]).norm()) #could cause error
        sol = ((data[i][2][0])/divider)
        Q_matrix = Q_matrix.col_insert(i, sol)
    lambda_matrix = diag(*egenvals)

    if text == False:
        return Q_matrix, lambda_matrix
    
    string1 = r'''Da \: A \: er \: symetrisk \: findes \: der \: en \: ortogonal \: matrix, \: Q, \: således \: at \: Q^{T}AQ=\Lambda '''
    string2 = r'''Lambda-matricen \: dannes \: ved \: at \: lade \: {%s}'s \: egenværdier \: løbe \: i \: diagonalen: \Lambda = {%s}''' % (retrieve_name(A)[0], latex(lambda_matrix))
    string3 = r'''Q-matricen, \: dannes \: ved \: at \: indsætte \: {%s}'s \: egenvektorer \: (normeret), \: som \: søjler \: i \: Q. \: Q = {%s}''' % (retrieve_name(A)[0], latex(Q_matrix))
    display(Math(string1))
    display(Math(string2))
    display(Math(string3))
    return Q_matrix, lambda_matrix


#Tager en funktion af 2 variable som input (kun x,y for nu)
def kvadratisk_form_til_reducerede_from(f, text=True):
    x1, y1 = symbols('x1,y1')
    x_x = f.coeff(x**2)
    y_y = f.coeff(y**2)
    x_y = f.coeff(x*y)
    kvadratisk_form_var = (x_x*x**2)+(y_y*y**2)+(x_y*x*y)
    A = Matrix([[x_x, x_y/2], [x_y/2, y_y]])
    Q_matrix, lambda_matrix = diagonalisering(A, text=False)
    reducerede_form = Matrix([x1,y1]).T*lambda_matrix*Matrix([x1,y1])
    reducerede_form_kun_det_kvadratiske = reducerede_form
    x_single = f.coeff(x).subs(y,0)
    y_single = f.coeff(y).subs(x,0)
    konstanter = f.subs(x,0).subs(y,0)
    reducerede_form += (Matrix([x_single, 0]).T*Q_matrix*Matrix([x1,y1]))
    reducerede_form += (Matrix([0, y_single]).T*Q_matrix*Matrix([x1,y1]))
    reducerede_form += Matrix([[konstanter]])

    if text == False:
        return reducerede_form
    
    string1 = r'''En \: del \: af \: funktionen, \: {%s}, \: optræder \: på \: kvadratisk \: form.\: Nemlig \: {%s}''' % (retrieve_name(f)[0], latex(kvadratisk_form_var))
    string2 = r'''Dette \: kan \: omskrives \: til: \: {%s}= \begin{matrix}x & y\end{matrix} A  \begin{matrix} x \\ y \end{matrix}, \: hvor \: A = {%s}''' % (latex(kvadratisk_form_var), latex(A))
    string3 = r'''Da \: A \: er \: symetrisk \: findes \: der \: en \: ortogonal \: matrix, \: Q, \: således \: at \: Q^{T}AQ=\Lambda '''
    string4 = r'''Gennem \: en \: række \: omregninger, \: kan \: det \: vises \: at \: \: \begin{matrix}x & y\end{matrix} A  \begin{matrix} x \\ y \end{matrix} = \begin{matrix} x_{1} \\ y_{1} \end{matrix}^{T}\Lambda \: \begin{matrix} x_{1} \\ y_{1} \end{matrix}={%s}''' % (latex(reducerede_form_kun_det_kvadratiske))
    string5 = r'''Den \: resterende \: del \: af \: funktionen, \: {%s}, \: kan \: ligeledes \: udtrykkes \: i \: det \: nye \: koordinatsysstem \: ved \: at \: anvende \: basisskiftematricen, \: Q''' % (retrieve_name(f)[0])
    string6 = r'''Hermed \: kan \: vi \: opskrive \: funktionen, \: {%s}, \: udtrykt \: i \: det \: nye \: koordinatsystem \: (x_{1},y_{1}) \: ved: {%s}''' % (retrieve_name(f)[0], latex(reducerede_form))
    string7 = r'''Basisskiftematrice, \: Q, \: er \: givet \: ved: Q = {%s}'''  % (latex(Q_matrix))
    display(Math(string1))
    display(Math(string2))
    display(Math(string3))
    display(Math(string4))
    display(Math(string5))
    display(Math(string6))
    display(Math(string7))
    
    return reducerede_form


# Tager et udtryk (udtrykt i x,y only!) og give complete square
# Et eksempel kunne være expr = 5*x**2 + 4*y**2 + 20*x + 12*y + 9
def complete_square(expr):
    x_x = expr.coeff(x**2)
    y_y = expr.coeff(y**2)
    x_ = expr.coeff(x)
    y_ = expr.coeff(y)
    complete_square_sol = x_x*(x+((x_/2)/x_x))**2 + y_y*(y+((y_/2)/y_y))**2
    complete_square_konstant = (expand(complete_square_sol)).subs(x,0).subs(y,0)
    expr_konstant = expr.subs(x,0).subs(y,0)

    if complete_square_konstant == expr_konstant:
        return complete_square_sol
    else:
        return Eq(complete_square_sol, abs(complete_square_konstant-expr_konstant))


#Bestemmer fluxen gennem et lukket rumligt område
def gauss(V, r, uint, vint, wint, text=True):
    integrant = divv(V, (x,y,z)).subs(x,r[0]).subs(y,r[1]).subs(z,r[2])
    sol = integrate(integrant*jacobi_rum(r), (u, uint[1], uint[2]), (v, vint[1], vint[2]), (w, wint[1], wint[2]))
    if text == False:
        return sol
    string1 = r'''Gauss' \: divergens-sætninger \: siger, \: at \: fluxen \int_{\partial \Omega}^{}V\cdot n \: d\mathit{\mu}, \: kan \: bestemmes \: ved \int_{\partial \Omega}^{}V\cdot n \: d\mathit{\mu} = \int_{\Omega r}^{}\mathit{Div}\left(V\right)d\mathit{\mu}, \: for \: et \: lukket \: rumligt \: område, \: \Omega'''
    string2 = r'''\int_{\partial \Omega}^{}V\cdot n \: d\mathit{\mu} = \int_{\Omega r}^{}\mathit{Div}\left(V\right)d\mathit{\mu} = \int_{{%s}}^{{%s}}\int_{{%s}}^{{%s}}\int_{{%s}}^{{%s}}\mathit{Div}\left(V\right)\left(r\left(u,v,w\right)\right) \mathit{Jacobi}\left(u,v,w\right)\mathit{dudvdw}''' % (latex(uint[1]), latex(uint[2]), latex(vint[1]), latex(vint[2]), latex(wint[1]), latex(wint[2]))
    string3 = r'''Divergensen \: er \: givet \: ved \: Div(V)={%s}''' % (latex(divv(V, (x,y,z))))
    string4 = r'''\mathit{Div}\left(V\right)\left(r\left(u,v,w\right)\right) = {%s}''' % (latex(integrant))
    string5 = r'''Jacobi_{r} = {%s}''' % (latex(jacobi_rum(r)))
    string6 = r'''Hermed \: er \: den \: samlede \: integrant \: givet \: ved: {%s}''' % (latex(integrant*jacobi_rum(r)))
    string7 = r'''Nu \: kan \: integrallet \: udregnes: \int_{{%s}}^{{%s}}\int_{{%s}}^{{%s}}\int_{{%s}}^{{%s}}\mathit{Div}\left(V\right)\left(r\left(u,v,w\right)\right) \mathit{Jacobi}\left(u,v,w\right)\mathit{dudvdw} = \int_{{%s}}^{{%s}}\int_{{%s}}^{{%s}}\int_{{%s}}^{{%s}} {%s} = {%s}''' % (latex(uint[1]), latex(uint[2]), latex(vint[1]), latex(vint[2]), latex(wint[1]), latex(wint[2]), latex(uint[1]), latex(uint[2]), latex(vint[1]), latex(vint[2]), latex(wint[1]), latex(wint[2]), latex(integrant*jacobi_rum(r)), latex(sol))
    display(Math(string1))
    display(Math(string2))
    display(Math(string3))
    display(Math(string4))
    display(Math(string5))
    display(Math(string6))
    display(Math(string7))
    return sol

def jacobi_plan(r, u, v, text=True):
    sol = abs(det(Matrix([[diff(r[0], u), diff(r[0], v)], [diff(r[1], u), diff(r[1], v)]])))
    if text == False:
        return sol
    string1 = r'''Jacobi-funktionen \: for \: et \: plant \: område \: i \: planen \: er \: givet \: ved: |\det\left(\left[\begin{array}{cc}
\frac{\partial x}{\partial u}r\left(u,v\right) & \frac{\partial x}{\partial v}r\left(u,v\right) 
\\
 \frac{\partial y}{\partial u}r\left(u,v\right) & \frac{\partial y}{\partial v}r\left(u,v\right) 
\end{array}\right]\right)| = |\det\left(\left[\begin{array}{cc}
{%s} & {%s} 
\\
 {%s} & {%s} 
\end{array}\right]\right)| = {%s}''' % (latex(diff(r[0], u)), latex(diff(r[0], v)), latex(diff(r[1], u)), latex(diff(r[1], v)), latex(sol))
    display(Math(string1))
    return sol


def planintegral_af_func_over_plan(f, r, uint, vint, text=True):
    fruv = f.subs(x, r[0]).subs(y, r[1])
    sol = integrate(fruv*jacobi_plan(r, u, v, text=False), (u, uint[1], uint[2]), (v, vint[1], vint[2]))
    if text == None:
        return simplify(sol)
    string1 = r'''Planintegralet \: af \: funktionen \: f\left(x,y\right) \: over \: det \: parametriserede \: område \: P_{r} \: defineres \: ved:'''
    string2 = r'''\int_{F_{r}} f \: d\mu = \int_{c}^{d}\int_{a}^{b}f\left(r\left(u,v\right)\right)\mathit{Jacobi}_{r}\left(u,v\right)\mathit{dudv}'''
    string3 = r'''Her \: er \: \mathit{Jacobi}_{r}\left(u,v\right) \: givet \: ved:'''
    string4 = r'''\mathit{Jacobi}_{\boldsymbol{\mathit{r}}}\left(u,v\right)=|\det\left(\left[\begin{array}{cc}
\frac{\partial x}{\partial u}r\left(u,v\right) & \frac{\partial x}{\partial v}r\left(u,v\right) 
\\
 \frac{\partial y}{\partial u}r\left(u,v\right) & \frac{\partial y}{\partial v}r\left(u,v\right) 
\end{array}\right]\right)| '''
    string5 = r'''Jeg \: bestemmer \: først \: f\left(r\left(u,v\right)\right) \: og \: \mathit{Jacobi}_{\boldsymbol{\mathit{r}}}\left(u,v\right) '''
    string6 = r'''f\left(r\left(u,v\right)\right) = {%s}''' % (latex(fruv))
    string7 = r'''\mathit{Jacobi}_{\boldsymbol{\mathit{r}}}\left(u,v\right)=| \: r_{u}^{\prime}\left(u,v\right)\times r_{v}^{\prime}\left(u,v\right)| = {%s} ''' % (latex(jacobi_plan(r, u, v, text=False)))
    string8 = r'''Nu \: indsættes \: i \: formlen \: sammen \: med \: grænserne \: u \in [%s ,%s]\mathit{og} \: v \in [%s ,%s]''' % (latex(uint[1]), latex(uint[2]), latex(vint[1]), latex(vint[2]))
    string9 = r'''\int_{P_{r}} f \: d\mu = \int_{%s}^{%s}\int_{%s}^{%s}{%s}{%s}\mathit{dudv} = {%s}''' % (latex(uint[1]), latex(uint[2]), latex(vint[1]), latex(vint[2]), latex(fruv), latex(jacobi_plan(r, u, v, text=False)), latex(sol))
    display(Math(string1))
    display(Math(string2))
    display(Math(string3))
    display(Math(string4))
    display(Math(string5))
    display(Math(string6))
    display(Math(string7))
    display(Math(string8))
    display(Math(string9))
    return simplify(sol)


def rumintegral_af_func_over_rumobjekt(f, r, uint, vint, wint, text=True):
    fruv = f.subs(x, r[0]).subs(y, r[1]).subs(z, r[2])
    sol = integrate(fruv*jacobi_rum(r, u, v, w, text=False), (u, uint[1], uint[2]), (v, vint[1], vint[2]), (w, wint[1], wint[2]))
    if text == None:
        return simplify(sol)
    string1 = r'''Rumintegralet \: af \: funktionen \: f\left(x,y,z\right) \: over \: det \: parametriserede \: rumlige \: område \: \Omega_{r} \: defineres \: ved:'''
    string2 = r'''\int_{\Omega_{r}} f \: d\mu = \int_{h}^{l}\int_{c}^{d}\int_{a}^{b}f\left(r\left(u,v,w\right)\right)\mathit{Jacobi}_{r}\left(u,v,w\right)\mathit{dudvdw}'''
    string3 = r'''Her \: er \: \mathit{Jacobi}_{r}\left(u,v,w\right) \: givet \: ved:'''
    string4 = r'''\mathit{Jacobi}_{\boldsymbol{\mathit{r}}}\left(u,v,w\right)=|\det\left(\left[\begin{array}{cc}
\frac{\partial x}{\partial u}r\left(u,v,w\right) & \frac{\partial x}{\partial v}r\left(u,v,w\right) & \frac{\partial x}{\partial w}r\left(u,v,w\right) 
\\
 \frac{\partial y}{\partial u}r\left(u,v,w\right) & \frac{\partial y}{\partial v}r\left(u,v,w\right) & \frac{\partial y}{\partial w}r\left(u,v,w\right)
\\
 \frac{\partial z}{\partial u}r\left(u,v,w\right) & \frac{\partial z}{\partial v}r\left(u,v,w\right) & \frac{\partial z}{\partial w}r\left(u,v,w\right)
\end{array}\right]\right)| '''
    string5 = r'''Jeg \: bestemmer \: først \: f\left(r\left(u,v,w\right)\right) \: og \: \mathit{Jacobi}_{\boldsymbol{\mathit{r}}}\left(u,v,w\right) '''
    string6 = r'''f\left(r\left(u,v,w\right)\right) = {%s}''' % (latex(fruv))
    string7 = r'''\mathit{Jacobi}_{\boldsymbol{\mathit{r}}}\left(u,v,w\right)=| \: r_{u}^{\prime}\left(u,v\right)\times r_{v}^{\prime}\left(u,v\right)| = {%s} ''' % (latex(jacobi_rum(r, u, v, w, text=False)))
    string8 = r'''Nu \: indsættes \: i \: formlen \: sammen \: med \: grænserne \: u \in [%s ,%s]\mathit{og} \: v \in [%s ,%s]\mathit{og} \: w \in [%s ,%s]''' % (latex(uint[1]), latex(uint[2]), latex(vint[1]), latex(vint[2]), latex(wint[1]), latex(wint[2]))
    string9 = r'''\int_{\Omega_{r}} f \: d\mu = \int_{%s}^{%s}\int_{%s}^{%s}\int_{%s}^{%s}{%s}{%s}\mathit{dudvdw} = {%s}''' % (latex(uint[1]), latex(uint[2]), latex(vint[1]), latex(vint[2]), latex(wint[1]), latex(wint[2]), latex(fruv), latex(jacobi_plan(r, u, v, text=False)), latex(sol))
    display(Math(string1))
    display(Math(string2))
    display(Math(string3))
    display(Math(string4))
    display(Math(string5))
    display(Math(string6))
    display(Math(string7))
    display(Math(string8))
    display(Math(string9))
    return simplify(sol)


def length_af_kurve(r, uint, text=True):
    sol = integrate(jacobi_kurve(r, u, text=False), (u, uint[1], uint[2]))
    if text==False:
        return sol
    string1 = r'''Længden \: af \: den \: parametriserede \: kurve \: K_{r} : r\left(u\right)=\left(x\left(u\right),y\left(u\right),z\left(u\right)\right)u\in \left[a,b\right]'''
    string2 = r'''defineres \: som \: kurveintegralet: \int_{K_{r}}^{}1d \mu = \int_{{a}}^{b}|r'\left(u\right)|d u'''
    string3 = r'''Hvor \: integranten \: netop \: er \: Jacobi \: funktionen \: \mathit{Jacobi}_{r}\left(u\right)=|r'\left(u\right)|'''
    string4 = r'''\mathit{Jacobi}_{r}\left(u\right)=|r'\left(u\right)|={%s}''' % (latex(jacobi_kurve(r, u, text=False)))
    string5 = r'''\int_{K_{r}}^{}1d \mu = \int_{{{%s}}}^{{%s}}{%s}d u = {%s}''' % (latex(uint[1]), latex(uint[2]), latex(jacobi_kurve(r, u, text=False)), latex(sol))
    display(Math(string1))
    display(Math(string2))
    display(Math(string3))
    display(Math(string4))
    display(Math(string5))
    return sol



def areal_af_plant_area(r, uint, vint, text=True):
    sol = integrate(jacobi_plan(r, u, v, text=False), (u, uint[1], uint[2]), (v, vint[1], vint[2]))
    if text==False:
        return sol
    string1 = r'''Arealet \: af \: et \: parametriseret \: plant \: område \: B_{r} : r\left(u,v\right)=\left(x\left(u,v\right),y\left(u,v\right)\right)u\in \left[a,b\right] v\in \left[c,d\right]'''
    string2 = r'''defineres \: som \: planintegralet: \int_{B_{r}}^{}1d \mu = \int_{{c}}^{d}\int_{{a}}^{b}Jacobi_{r} \:d u d v'''
    string4 = r'''\int_{B_{r}}^{}1d \mu = \int_{{{%s}}}^{{%s}}\int_{{{%s}}}^{{%s}}{%s}d u d v = {%s}''' % (latex(vint[1]), latex(vint[2]), latex(uint[1]), latex(uint[2]), latex(jacobi_plan(r, u, v, text=False)), latex(sol))
    display(Math(string1))
    display(Math(string2))
    string3 = r'''Hermed \: er \: Jacobi \: funktionen \: givet \: ved: {%s}''' % (latex(jacobi_plan(r, u, v)))
    display(Math(string3))
    display(Math(string4))
    return sol


def fladeintegral_areal(r, uint, vint, text=True):
    sol = integrate(jacobi_flade(r, u ,v, text=False), (u, uint[1], uint[2]), (v, vint[1], vint[2]))
    if text==False:
        return sol
    string1 = r'''Arealet \: af \: den \: parametriserede \: flade \: F_{r} : r\left(u,v\right)=\left(x\left(u,v\right),y\left(u,v\right),z\left(u,v\right)\right)u\in \left[a,b\right] v\in \left[c,d\right]'''
    string2 = r'''defineres \: som \: fladeintegralet: \int_{F_{r}}^{}1d \mu = \int_{{c}}^{d}\int_{{a}}^{b}Jacobi_{r} \:d u d v'''
    display(Math(string1))
    display(Math(string2))
    string3 = r'''Hermed \: er \: Jacobi \: funktionen \: givet \: ved: {%s}''' % (latex(jacobi_flade(r, u, v)))
    display(Math(string3))
    string4 = r'''\int_{F_{r}}^{}1d \mu = \int_{{{%s}}}^{{%s}}\int_{{{%s}}}^{{%s}}{%s}d u d v = {%s}''' % (latex(vint[1]), latex(vint[2]), latex(uint[1]), latex(uint[2]), latex(jacobi_flade(r, u, v, text=False)), latex(sol))
    display(Math(string4))

    return sol



def kurveintegral_af_f_over_kurve(f, r, u, start, slut, text=True):
    fru = f.subs(x, r[0]).subs(y, r[1]).subs(z, r[2])
    sol = integrate(fru*jacobi_kurve(r, u, text=False), (u, start, slut))
    if text == False:
        return sol
    string1 = r'''Kurveintegralet \: af \: funktionen \: f\left(x,y,z\right) \: over \: det \: parametriserede \: kurve \: K_{r} \: defineres \: ved:'''
    string2 = r'''\int_{K_{r}} f \: d\mu = \int_{a}^{b}f\left(r\left(u\right)\right)\mathit{Jacobi}_{r}\left(u\right)\mathit{du}'''
    string4 = r'''Jeg \: bestemmer \: nu \: f\left(r\left(u\right)\right):'''
    string5 = r'''f\left(r\left(u\right)\right) = {%s}''' % (latex(fru))
    string6 = r'''Nu \: indsættes \: i \: formlen \: sammen \: med \: grænserne \: u \in [%s ,%s]''' % (latex(start), latex(slut))
    string7 = r'''\int_{K_{r}} f \: d\mu = \int_{%s}^{%s}{%s}{%s}\mathit{du} = {%s}''' % (latex(start), latex(slut), latex(fru), latex(jacobi_kurve(r, u, text=False)), latex(sol))
    display(Math(string1))
    display(Math(string2))
    string3 = r'''Her \: er \: \mathit{Jacobi}_{r}\left(u\right) \: altså \: givet \: ved: {%s}''' % (latex(jacobi_kurve(r,u)))
    display(Math(string3))
    display(Math(string4))
    display(Math(string5))
    display(Math(string6))
    display(Math(string7))
    return sol