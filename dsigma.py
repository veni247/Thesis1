import numpy as np
import matplotlib.pyplot as plt

# Định nghĩa các hàm
def gV(e, sw, cw):
    return (e / (2 * sw * cw)) * (4 * sw**2 - 1)

def gA(e, sw, cw):
    return - e / (2 * sw * cw)

def s(E):
    return 16 * E**4

def abs_chi0_squared(s, e, m_Z, Gamma_Z):
    return s**2 / (16 * e**4 * ((s - m_Z**2)**2 + Gamma_Z**2 * m_Z**2))

# Tính phần thực của χ₀
def real_chi0(s, e, m_Z, Gamma_Z):
    return (s - m_Z**2) / ((s - m_Z**2)**2 + Gamma_Z**2 * m_Z**2)

# Định nghĩa các hàm G₁(s), G₂(s), và G₃(s)
def G1(s, gV, gA, e, m_Z, Gamma_Z):
    abs_chi0_sq = abs_chi0_squared(s, e, m_Z, Gamma_Z)
    real_chi0_val = real_chi0(s, e, m_Z, Gamma_Z)
    return abs_chi0_sq * (gV**2 + gA**2)**2 + 2 * gV**2 * real_chi0_val + 1

def G2(s, gV, gA, e, m_Z, Gamma_Z):
    abs_chi0_sq = abs_chi0_squared(s, e, m_Z, Gamma_Z)
    real_chi0_val = real_chi0(s, e, m_Z, Gamma_Z)
    return 2 * gA**2 * gV**2 * abs_chi0_sq + gA**2 * real_chi0_val

def G3(s, gV, gA, e, m_Z, Gamma_Z, mu):
    abs_chi0_sq = abs_chi0_squared(s, e, m_Z, Gamma_Z)
    real_chi0_val = real_chi0(s, e, m_Z, Gamma_Z)
    return (abs_chi0_sq * ((gV**2 + gA**2)**2 + (gV**4 - gA**4) * mu) 
            + 2 * gV**2 * real_chi0_val * (mu + 1) + mu + 1)

# Hàm vi phân mật độ phân tán dσ/dΩ
def dsigma_dOmega(e, k, E, s, theta, G1_val, G2_val, G3_val):
    cos_theta = np.cos(theta)
    return np.sin(theta) * (e**4 * abs(k) / (8 * np.pi * E * s**2) 
            * (G1_val * abs(k)**2 * cos_theta**2 + 4 * E * G2_val * abs(k) * cos_theta + G3_val * E**2))

# Hằng số
e = 1.6  # Điện tích cơ bản (Coulombs)
sw = 0.23  # sin²(theta_W)
cw = np.sqrt(1 - sw**2)  # cos(theta_W)
m_Z = 90  # Khối lượng boson Z (GeV)
Gamma_Z = 2.4952  # Độ rộng của boson Z (GeV)
m_u = 0.105  # Khối lượng quark lên (GeV)

# Đặt giá trị √s
sqrt_s_values = [200, 20]
s_values = np.array(sqrt_s_values) ** 2

# Tạo hình và trục
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Vòng lặp qua từng giá trị sqrt(s)
for idx, sqrt_s in enumerate(sqrt_s_values):
    s_val = sqrt_s ** 2
    E = (s_val / 16)**0.25  # Giải cho E

    # Tính gV, gA và giá trị G
    gV_val = gV(e, sw, cw)
    gA_val = gA(e, sw, cw)

    # Tính G1, G2, G3
    G1_val = G1(s_val, gV_val, gA_val, e, m_Z, Gamma_Z)
    G2_val = G2(s_val, gV_val, gA_val, e, m_Z, Gamma_Z)
    G3_val = G3(s_val, gV_val, gA_val, e, m_Z, Gamma_Z, m_u**2 / E**2)

    # Động lượng của hạt phát ra k
    k = np.sqrt(E**2 - m_u**2)

    # Phạm vi các giá trị theta
    theta_values = np.linspace(0, np.pi, 500)

    # Tính mật độ phân tán
    dsigma_values = [dsigma_dOmega(e, k, E, s_val, theta, G1_val, G2_val, G3_val) for theta in theta_values]

    # Vẽ kết quả cho từng trục
    axes[idx].plot(theta_values, dsigma_values, label=f"$\\sqrt{{s}} = {sqrt_s} \\, GeV$")
    axes[idx].set_xlabel(r"$\theta$ (radian)")
    axes[idx].set_ylabel(r"$\frac{d\sigma}{d\theta}$")
    axes[idx].legend()
    axes[idx].grid()
    axes[idx].set_title(f"$\\sqrt{{s}} = {sqrt_s} \\, GeV$")

# Hoàn thiện biểu đồ
plt.suptitle("Mật độ phân tán vi phân $d\\sigma/d\\theta$ theo góc tán xạ $\\theta$")
plt.tight_layout(rect=[0, 0, 1, 0.96])  # Điều chỉnh để có không gian cho tiêu đề
plt.show()
