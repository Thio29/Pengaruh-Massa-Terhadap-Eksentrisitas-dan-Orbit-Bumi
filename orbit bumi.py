# Import library
import numpy as np
import matplotlib.pyplot as plt

# Fungsi Program
def program():
    print()
    print("Program Analisa Pengaruh Massa Matahari Terhadap Bentuk".center(60, " "))
    print("dan Eksentrisitas Orbit Bumi Mengelilingi Matahari".center(60, " "))
    menu()
    
# Fungsi Menu
def menu():
    
    # Fungsi jumlah data
    def jumlah_data(menu):
            
        # penentuan minimal banyak data
        if menu == 1:
            minim = 1
        elif menu == 2:
            minim = 2
            
        try:
            n = int(input("\nJumlah Data : "))
        except ValueError:
            print("\nHarap masukkan bilangan bulat. Silahkan ulangi")
            return jumlah_data(menu)
            
        if n < minim or n > 10:
            print("\nHarap masukkan jumlah data sebanyak %d─10. Silahkan ulangi."%minim)
            return jumlah_data(menu)
            
        return n
    
    # Fungsi input massa
    def input_massa(n):
        
        # List massa    
        M = []
            
        print("\nInput massa matahari (default M = 2x10³⁰kg):")
        
        # Iterasi untuk menginput data massa    
        for i in range(1, n+1):
            try:
                massa = float(input("Massa-%d (kg): " %i))
            except ValueError:
                print("\nHarap masukkan angka (bukan teks). Silahkan ulangi")
                return input_massa(n)
            if massa < 2e30 or massa > 20e30:
                print("\nHarap masukkan besar massa matahari pada rentang \n2x10³⁰ ≤ M ≤ 20x10³⁰ kg. Silahkan ulangi.")
                return input_massa(n)
            else:
                M.append(massa)
                
        return np.sort(np.array(M))
        
    # Fungsi menu 1    
    def menu1():
        n = jumlah_data(1)
        M = input_massa(n)
        bentuk_orbit(M)
        
    # Fungsi menu 2
    def menu2():
        n = jumlah_data(2)
        M = input_massa(n)
        e = eksentrisitas(M)
        M = M / 1e30
        regresi(M, e)
    
    # Fungsi menu akhir    
    def menuakhir():
        print()
        print("1. Kembali ke menu")
        print("2. Akhiri program")
        try:
            pilih = int(input('\nPilih Menu (1/2):'))
        except ValueError:
            print("\nHarap masukkan bilangan bulat. Silahkan ulangi.")
            menuakhir()
            quit()
        if pilih == 1:
            menu()
        elif pilih == 2:
            keluar()
        else:
            print("\nHarap masukkan pilihan menu yang ada. Silahkan ulangi")
            menuakhir()
    
    # Fungsi keluar        
    def keluar():
        print()
        print("Terima kasih telah menggunakan program kami!".center(60, " "))
        print("Sampai jumpa".center(60, "-"))
        exit()
    
    # Menu awal
    print()
    print(" Menu Program ".center(60, "-"))
    print("1.  Simulasi Pengaruh Perubahan Massa Matahari Terhadap\n    Bentuk Orbit Bumi")
    print("2.  Analisa Pengaruh Perubahan Massa Matahari Terhadap\n    Eksentrisitas Orbit Bumi")
    print("99. Akhiri program")
    try:
        pilih = int(input('\nPilih Menu (1/2/99):'))
    except ValueError:
        print("\nHarap masukan bilangan bulat. Silahkan Ulangi")
        menu()
    if pilih == 1:
        menu1()
        menuakhir()
    elif pilih == 2:
        menu2()
        menuakhir()
    elif pilih == 99:
        keluar()
    else:
        print("\nHarap masukan pilihan menu yang ada. Silahkan Ulangi")
        menu()

# Penyelesaian PDB
def masalah_nilai_awal(M):
        
    # Periode Bumi Mengelilingi Matahari
    def periode(x0, M):
        G = 6.67408e-11             # konstanta umum gravitasi m^3/kg s^2
        r = np.sqrt(x0[0]**2+x0[1]**2+x0[2]**2)
        T = np.sqrt((4 * np.pi**2 * r**3)/(G*M))
        return T
        
    # Persamaan Gerak Bumi
    def modelbumi(y, t, m):
        G = 6.67408e-11             # konstanta umum gravitasi m^3/kg s^2   
        x = y[0:3].copy()
        v = y[3:6].copy()
        r = np.sqrt(x[0]**2+x[1]**2+x[2]**2)
        dxdt = v
        dvdt = -G * m * x / r**3
        dy = np.hstack((dxdt,dvdt))
        return dy
        
    # Metode Runge-Kutta orde 4
    def rungekutta4(f, y0, t, m):
        Y = []
        Y.append(y0)
        h = t[1]-t[0]
        for i in range(len(t)):
            k1 = h * (f(y0, t[i], m))
            k2 = h * (f((y0+k1/2), (t[i]+h/2), m))
            k3 = h * (f((y0+k2/2), (t[i]+h/2), m))
            k4 = h * (f((y0+k3), (t[i]+h), m))
            k = (k1+2*k2+2*k3+k4)/6
            yn = y0 + k
            y0 = yn
            Y.append(yn)
        return np.array(Y)
    
    # Nilai astronomi
    au = 1.496e11                   # satuan astronomi [m]
    
    # Posisi awal
    x0 = au
    r0 = np.array([x0, 0., 0.])
    
    # Kecepatan awal
    vx0 = 29290 
    v0 = np.array([0., vx0, 0.])
    
    # Keadaan awal
    s0 = np.hstack((r0, v0))
      
    # Waktu awal
    t0 = 0.0
    
    # Banyak Step
    nt = 10000
        
    # iterasi penyelesaian PDB
    T = periode(r0, M)
    t = np.linspace(t0, T, nt)
    sol = rungekutta4(modelbumi, s0, t, M)
    x = sol[:,0:3].T
    v = sol[:,3:6].T
    return x, v

# Bentuk orbit    
def bentuk_orbit(M):
    # Nilai astronomi
    au = 1.496e11                   # satuan astronomi [m]
    
    # Mengatur grafik
    plt.figure()
    
    # Mendapatkan solusi PDB
    for i in range(len(M)):
        x, v = masalah_nilai_awal(M[i])
        m = np.format_float_scientific(M[i])
        
        # Plotting orbit
        if i == 0:
            plt.plot(x[0,:]/au, x[1,:]/au,label="M = %s kg"%m)
        else:
            plt.plot(x[0,:]/au, x[1,:]/au,linestyle=":",label="M = %s kg"%m)
    
    # Pengaturan grafik
    plt.plot([0.],[0.],color='orange', marker="o")
    plt.xlabel('x (satuan astronomi)')
    plt.ylabel('y (satuan astronomi)')
    plt.title("Bentuk Orbit Bumi mengelilingi Matahari (2D)")
    plt.axis('equal')
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.show()
    
# Eksentrisitas orbit   
def eksentrisitas(M):
    
    # list eksentrisitas
    ek = []
    
    # Mendapatkan solusi PDB
    for i in range(len(M)):
        x, v = masalah_nilai_awal(M[i])
        
        # mendapatkan nilai radius
        ra = abs(max(x[0]))  # apoapsis
        rp = abs(min(x[0]))  # periapsis
        
        # mendapatkan nilai eksentrisitas
        e = 1 - 2/((ra/rp)+1)
        
        # memasukkan nilai e ke dalam list
        ek.append(e)
    
    return np.array(ek)
  
# Metode Regresi
def regresi(x, y):
    
    # Standar deviasi
    def standar_deviasi(m, y, y1):
        n = len(y)
        sigma = 0.0
        for i in range(n):
            sigma += (y[i] - y1[i])**2
        std = np.sqrt(sigma/(n-m))
        return std  
    
    # Regresi polinomial
    def polinomial(x, y, m):
        
        # Menentukan persamaan normal polinomial
        def polyfit(x,y,m):
            a = np.zeros((m+1, m+1))
            b = np.zeros(m+1)
            s = np.zeros(2*m+1)
            for i in range(len(x)):
                temp = y[i]
                for j in range(m+1):
                    b[j] = b[j] + temp
                    temp = temp*x[i]
                temp = 1.0
                for j in range(2*m+1):
                    s[j] = s[j] + temp
                    temp = temp*x[i]
            for i in range(m+1):
                for j in range(m+1):
                    a[i,j] = s[i+j]
            return a, b
        
        # Menentukan solusi persamaan normal
        def gauss_seidel(x, A, B):
            n = len(x)
            E = 1
            while E > 1e-6:
                xb = []
                for i in range(n):
                    sigma = 0
                    for j in range(0,i):
                        sigma = sigma + A[i][j]*xb[j]
                    for j in range(i+1,n):
                        sigma = sigma + A[i][j]*x[j]
                    hasil = (B[i]-sigma)/A[i][i]
                    xb.append(hasil)
                E = abs(max(x)-max(xb))
                x = np.array(xb)
            return x
        
        # Mendapatkan koefisien polinomial
        A, B = polyfit(x, y, m)
        z = np.zeros_like(B)
        c = gauss_seidel(z, A, B)
        return c
    
    # Regresi linier
    def linier(x, y):
        # Rata-rata x dan y
        xr = sum(x) / len(x)
        yr = sum(y) / len(y)
        # menghitung konstanta b
        num = 0; denom = 0
        for i in range(len(x)):
            num = num + y[i]*(x[i]-xr)
            denom = denom + x[i]*(x[i]-xr)
        b = num/denom
        # menghitung koefisien a
        a = yr - xr*b
        return a, b
        
    # dictionary regresi
    regresi = {}
    
    # regresi linier
    a, b = linier(x, y)
    yr = np.array(a+b*x)
    stdr = standar_deviasi(1, y, yr)
    regresi[0] = {"nama": "linier", "koefisien" : np.array([a, b]), "hampiran": yr, "std": stdr}
    
    # menentukan batas akhir iterasi
    if len(x) > 5:
        n = 4
    elif len(x) <= 5:
        n = len(x)
    
    # mendapatkan koefisien
    for m in range(1, n):
        c = polinomial(x, y, m)
        
        # mendapatkan nilai hampiran
        yp =np.zeros(len(x))*1.0
        for i in range(m+1):
            yp = yp + c[i]*x**i
        
        # standar deviasi
        std = standar_deviasi(m, y, yp)
        
        # menambahkan nilai ke ditionary
        regresi[m] = {"nama": "polinomial derajat %d"%m, "koefisien" : c, "hampiran": yp, "std": std}
    
    # mencari standar deviasi terkecil dari tiap persamaan
    std = []
    for i in range(len(regresi)):
        std.append(regresi[i]["std"])
    
    # mendapatkan persamaan terbaik
    for i in range(len(std)):
        if min(std) == regresi[i]["std"]:
            c = np.array(regresi[i]["koefisien"])
            
            # mendapatkan persamaan
            if i == 0 or i == 1:
                eq = "%fx + %f" %(c[1], c[0])
            elif i == 2:
                eq = "%fx² + %fx + %f" %(c[2], c[1], c[0])
            elif i == 3:
                eq = "%fx³ + %fx² + %fx + %f" %(c[3], c[2], c[1], c[0])
            elif i == 4:
                eq = "%fx⁴ + %fx³ + %fx² + %fx + %f" %(c[4], c[3], c[2], c[1], c[0])
      
    # mendapatkan titik data kurva (x dan y)
    xplot = np.linspace(min(x), max(x), 1000)
    yplot = np.zeros(len(xplot))*1.0
    for i in range(len(c)):
        yplot = yplot + c[i]*xplot**i
    
    # plotting data dan kurva
    plt.figure()
    plt.plot(x, y, "o")
    plt.plot(xplot, yplot, label=eq)
    plt.xlabel("Massa matahari (10³⁰ kg)")
    plt.ylabel("Eksentrisitas Orbit")
    plt.title("Hubungan Massa Matahari dengan Eksentrisitas Orbit Bumi")
    plt.legend(loc=2)
    plt.grid()
    plt.tight_layout()
    plt.show()
    
# Program Utama
if __name__ == '__main__':
    program()
    
