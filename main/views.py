from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib.auth.views import LoginView
from main.static.tubes_calculate import *
from .forms import *
from .static.tubes_calculate import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from io import StringIO


x = [[],[],[],[]]
y = [[],[],[],[]]
z = [[],[],[],[]]
lx1 = []
lx2 = []
ly1 = []
ly2 = []
lz1 = []
lz2 = []
angle = 45
walls = []

def return_graph():
    global x, y, z
    global lx1, lx2, ly1, ly2, lz1, lz2
    global angle

    x_res = [[],[],[],[]]
    y_res = [[],[],[],[]]
    z_res = [[],[],[],[]]
    for index in range(len(x)):
        x_res[index] = ([i for i in x[index] if i is not None])
        y_res[index] = ([i for i in y[index] if i is not None])
        z_res[index] = ([i for i in z[index] if i is not None])

    # x = np.arange(0,np.pi*3,.1)
    # y = np.sin(x)

    fig = plt.figure(figsize=(7, 6))
    ax = plt.axes(projection='3d')
    #for i in range(len(x)):

    for i in range(len(lx1)):
        if lx1[i] is not None:
            ax.plot([lx1[i], lx2[i]], [ly1[i], ly2[i]], [lz1[i], lz2[i]], 'b', linewidth='4')

    r3, = ax.plot([], [], [], 'm', linewidth='5', label='Fitting')

    if len(x) != 0:
        if len(x[0]) != 0:
            for i in range(len(x[0])):
                if x[2][i] is None:
                    r3, = ax.plot([x[1][i], x[0][i]], [y[1][i], y[0][i]], [z[1][i], z[0][i]], '#FC4C02', linewidth='5', label='Fitting')
                else:
                    xc = (x[1][i] + x[0][i]) / 2
                    yc = (y[1][i] + y[0][i]) / 2
                    zc = (z[1][i] + z[0][i]) / 2

                    r3, = ax.plot([x[1][i], x[0][i]], [y[1][i], y[0][i]], [z[1][i], z[0][i]],  '#FC4C02' , linewidth='5', label='Fitting')
                    r3, = ax.plot([xc, x[2][i]], [yc, y[2][i]], [zc, z[2][i]], '#FC4C02', linewidth='5', label='Fitting')

                    if x[3][i] is not None:
                        r3, = ax.plot([xc, x[3][i]], [yc, y[3][i]], [zc, z[3][i]], '#FC4C02', linewidth='5', label='Fitting')

    r1, = ax.plot(x[0], y[0], z[0],'og', label='Fitting in port')
    r2, = ax.plot([],[],[],'or', label='Fitting out ports')

    for i in range(1, len(x_res)):
        r2, = ax.plot(x_res[i], y_res[i], z_res[i],'oc', label='Fitting out ports')




    plt.legend(handles = [r3, r1, r2])
    ax.view_init(30, angle)

    mpl.rcParams.update({"axes.grid": True, "grid.color": "#FFFFFF"})
    imgdata = StringIO()
    ax.set_xlabel('x[m]')
    ax.set_ylabel('y[m]')
    ax.axis('equal')

    plt.locator_params(nbins=5)

    fig.savefig(imgdata, format='svg')
    imgdata.seek(0)


    data = imgdata.getvalue()
    return data


def main(request):

    if request.method == "POST":
        form1 = Form1(request.POST)
        form2 = Form2(request.POST)
        form3 = Form3(request.POST)
        form4 = Form4(request.POST)
        form5 = Form5(request.POST)
        form6 = Form6(request.POST)
        form7 = Form7(request.POST)
        form8 = Form8(request.POST)
        form9 = Form9(request.POST)
        helpForm = HelpfulForm(request.POST)

        a = request.POST.get('GetMe', False)
        global angle
        if a:
            angle += float(a)
            return_graph()
        else:
            PU = []
            forms = [form1, form2, form3]
            for i in range(0, 3):
                form = forms[i]
                if form.is_valid():
                    if None not in form.cleaned_data.values():
                        PU.append((float(form.cleaned_data['X'+str(i+1)]), float(form.cleaned_data['Y'+str(i+1)]),
                                   float(form.cleaned_data['Z'+str(i+1)]), float(form.cleaned_data['D'+str(i+1)])))

            walls_forms = [form4, form5, form6]
            global walls
            walls = []
            Ds = []
            if form7.is_valid():
                for i in form7.cleaned_data.values():
                    if i is not None:
                        Ds.append(i)

            for i in range(len(walls_forms)):
                if walls_forms[i].is_valid() and len(Ds)>i:
                    if None not in walls_forms[i].cleaned_data.values():
                        walls.append((float(walls_forms[i].cleaned_data["X1"+str(i+1)]), float(walls_forms[i].cleaned_data["X2"+str(i+1)]),
                                      float(walls_forms[i].cleaned_data["Y1"+str(i+1)]), float(walls_forms[i].cleaned_data["Y2"+str(i+1)]), float(Ds[i])))

            # walls = [(-1.15, 1.568, -0.1535, -0.1535, 0.353),(1.568, 1.568, -0.1535, -2.710, 0.1), (1.568, -1.15, -2.710, -2.710, 0.05), (-1.15, -1.15, -0.1535, -1, 0.1)]
            # PU = [(0.34745, -0.153, 0.2485, 110), (1.567, -1.2301, 0.616, 50), (1.567, -1.943, 0.25850, 50), (1.585069, -2.648512, 0.410237, 40), (-1.14, -0.9, 0.610237, 50)]

            #PU = [(0.34745, -0.153, 0.2485, 110), (1.567, -1.2301, 0.616, 50), (1.567, -1.943, 0.4850, 50)]
            # boss = PlumbingUnit(1.585069, -2.648512, 0.410237, 40)

            # walls = [(-1.15, 1.568, -0.1535, -0.1535, 0.353),(1.568, 1.568, -0.1535, -2.710, 0.1), (1.568, -1.15, -2.710, -2.710, 0.05), (-1.15, -1.15, -0.1535, -1, 0.1)]
            # PU = [(0.34745, -0.153, 0.2485, 110), (1.567, -1.2301, 0.616, 50), (1.567, -1.943, 0.75850, 50)]

            d1 = {}
            d2 = {}
            d3 = {}
            room = initialize_room(PU, walls)

            import_tables('db.sqlite3', room)

            room.calculate_tubes()

            d1, d2, d3 = room.export_data('dataframe.xlsx')

            global x, y, z
            global lx1, lx2, ly1, ly2, lz1, lz2

            x = [d3['Вход. X'], d3['Выход0. X'], d3['Выход1. X'], d3['Выход2. X']]

            y = [d3['Вход. Y'], d3['Выход0. Y'], d3['Выход1. Y'], d3['Выход2. Y']]

            z = [d3['Вход. Z'], d3['Выход0. Z'], d3['Выход1. Z'], d3['Выход2. Z']]
            print(lx1)
            lx1 = d2['X1']
            lx2 = d2['X2']
            ly1 = d2['Y1']
            ly2 = d2['Y2']
            lz1 = d2['Z1']
            lz2 = d2['Z2']


            return HttpResponseRedirect('')

    else:
        form1 = Form1(request.POST)
        form2 = Form2(request.POST)
        form3 = Form3(request.POST)
        form4 = Form4(request.POST)
        form5 = Form5(request.POST)
        form6 = Form6(request.POST)
        form7 = Form7(request.POST)
        form8 = Form8(request.POST)
        form9 = Form9(request.POST)
        helpForm = HelpfulForm(request.POST)

    # rendered_form1 = form1.render("form_snippet.html")
    # rendered_form2 = form2.render("form_snippet.html")
    # rendered_form3 = form3.render("form_snippet.html")
    return render(request, 'index.html',
                  {'form1': form1, 'form2': form2, 'form3': form3, 'form4': form4, 'form5': form5, 'form6': form6,
                   'form7': form7, 'form8': form8, 'form9': form9, 'helpForm': helpForm, 'graph': return_graph()})
