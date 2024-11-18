set terminal pngcairo size 800,600 enhanced font 'Arial,10'
set output 'sigma_total.png'

set logscale y  # Đặt thang logarit cho trục y
set xlabel 'sqrt(s) (GeV)'
set ylabel 'sigma_total (pb)'
set title 'Đồ thị sigma_total theo sqrt(s)'
plot "sigma_total.txt" using 1:2 with lines title "$\sigma$ theo $\sqrt{s}$"