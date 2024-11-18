import numpy as np
import matplotlib.pyplot as plt

def gV(e, sw, cw):
    return (e / (2 * sw * cw)) * (4 * sw**2 - 1)

def gA(e, sw, cw):
    return - e / (2 * sw * cw)

def s(E):
    return 16 * E**4

def sigma_total(k, E, g_V, g_A, m_Z, Gamma_Z, m_mu, e):
    s_value = s(E)  # Tính giá trị s từ E

    # Tính toán các thành phần cần thiết
    term1 = (g_V**2 + g_A**2)**2 * s_value**2 / (16 * e**4 * ((s_value - m_Z**2)**2 + (Gamma_Z**2 * m_Z**2)))
    term2 = 2 * g_V**2 * s_value * (s_value - m_Z**2) / (4 * e**2 * ((s_value - m_Z**2)**2 + (Gamma_Z**2 * m_Z**2)))
    term3 = 1

    part1 = (term1 + term2 + term3) * (abs(k)**2) / 3

    term4 = (s_value**2 / (16 * e**4 * ((s_value - m_Z**2)**2 + (Gamma_Z**2 * m_Z**2)))) * (
        (g_V**2 + g_A**2)**2 + (g_V**4 - g_A**4) * (m_mu**2 / E**2)
    )
    term5 = (2 * g_V**2 * s_value * (s_value - m_Z**2) / (4 * e**2 * ((s_value - m_Z**2)**2 + (Gamma_Z**2 * m_Z**2)))) * (
        (m_mu**2 / E**2) + 1
    )
    term6 = (m_mu**2 / E**2) + 1

    part2 = (term4 + term5 + term6) * (E**2)

    # Tính toán sigma_total
    sigma_total_value = (abs(k) * e**4 / (4 * np.pi * E * s_value**2)) * (part1 + part2)
    
    return sigma_total_value

def main():
    # Tham số đầu vào
    sw = 0.23  # sin^2 theta_w
    cw = np.sqrt(1 - sw**2)  # cos^2 theta_w

    m_Z = 90  # Khối lượng boson Z (GeV)
    Gamma_Z = 2.4952  # Độ rộng của boson Z (GeV)
    m_mu = 0.105  # Khối lượng muon (GeV)

    # Tạo một dải E để tính toán
    E_values = np.linspace(0.1, 10, 500)  # Tăng số điểm để vẽ mượt hơn
    e = 1.6  # Hằng số điện tích của electron (Coulombs)

    # Tính toán sigma cho mỗi giá trị E
    sigma_values = []
    sqrt_s_values = []
    for E in E_values:
        g_V_value = gV(e, sw, cw)
        g_A_value = gA(e, sw, cw)
        k = np.sqrt(E**2 - m_mu**2)  # Đặt k = sqrt(E^2 - m_mu^2)
        sigma_value = sigma_total(k, E, g_V_value, g_A_value, m_Z, Gamma_Z, m_mu, e)
        sigma_values.append(sigma_value)
        sqrt_s_values.append(np.sqrt(s(E)))  # Tính giá trị sqrt(s)

    # Vẽ kết quả
    plt.figure(figsize=(10, 6))
    plt.plot(sqrt_s_values, sigma_values, label=r'$\sigma_{total}$', color='blue')
    plt.xlabel(r'$\sqrt{s}$ (GeV)')
    plt.xlim(0, 200)
    plt.xticks(np.arange(0, 201, 10))
    plt.ylabel(r'$\sigma_{total}$ (pb)')
    plt.title("Tiết diện tán xạ toàn phần")
    plt.grid()
    plt.legend()
    plt.yscale('log')  # Sử dụng thang logarit cho trục y
    plt.show()

    # Ghi kết quả vào file sigma_total.txt
    with open('sigma_total.txt', 'w', encoding='utf8') as f:
        # Ghi dữ liệu theo định dạng bảng, mỗi cột cách nhau bởi dấu cách
        for sqrt_s, sigma in zip(sqrt_s_values, sigma_values):
            f.write(f"{sqrt_s:.5f} {sigma:.5e}\n")  # Sử dụng dấu cách làm dấu phân cách

if __name__ == "__main__":
    main()
