# model.py - Базовая модель фильтрации нефти методом апвинд
import numpy as np


class OilFiltrationModel:
    """
    Базовая модель одномерной фильтрации нефти в пористой среде
    с использованием метода апвинд
    """

    def __init__(self):
        # Параметры пласта
        self.length = 100.0  # длина пласта, м
        self.porosity = 0.2  # пористость

        # Параметры флюидов
        self.mu_oil = 5.0  # вязкость нефти, мПа·с
        self.mu_water = 1.0  # вязкость воды, мПа·с
        self.initial_water_saturation = 0.2  # начальная водонасыщенность
        self.residual_oil_saturation = 0.2  # остаточная нефтенасыщенность

        # Параметры расчёта
        self.nx = 100  # число узлов сетки
        self.dx = self.length / self.nx  # шаг по x
        self.days = 100  # дней симуляции
        self.dt = 0.05  # шаг по времени, дней
        self.nt = int(self.days / self.dt) + 1  # число временных шагов

        # Создаем сетки
        self.x = np.linspace(0, self.length, self.nx + 1)
        self.t = np.linspace(0, self.days, self.nt)

        # Создаем массивы для хранения результатов
        # Насыщенность с учетом и без учета капиллярных эффектов
        self.Sw_with_cap = np.ones((self.nt, self.nx + 1)) * self.initial_water_saturation
        self.Sw_without_cap = np.ones((self.nt, self.nx + 1)) * self.initial_water_saturation

        # Устанавливаем граничные условия - закачка воды на входе
        self.Sw_with_cap[:, 0] = 0.8
        self.Sw_without_cap[:, 0] = 0.8

        # Параметры для модели капиллярного давления Брукса-Кори
        self.entry_pressure = 1.0  # давление входа, МПа
        self.pore_distribution_index = 1.5  # индекс распределения пор (λ)
        self.wettability_factor = 0.6  # коэффициент смачиваемости (1 - гидрофильная, 0 - гидрофобная)

    def relative_permeability_water(self, Sw):
        """Относительная проницаемость для воды"""
        Swc = self.initial_water_saturation
        Sor = self.residual_oil_saturation

        if Sw <= Swc:
            return 0.0
        elif Sw >= 1 - Sor:
            return 1.0
        else:
            Swn = (Sw - Swc) / (1 - Swc - Sor)
            return Swn ** 3  # кубическая зависимость

    def relative_permeability_oil(self, Sw):
        """Относительная проницаемость для нефти"""
        Swc = self.initial_water_saturation
        Sor = self.residual_oil_saturation

        if Sw >= 1 - Sor:
            return 0.0
        elif Sw <= Swc:
            return 1.0
        else:
            Son = (1 - Sw - Sor) / (1 - Swc - Sor)
            return Son ** 2  # квадратичная зависимость

    def fractional_flow(self, Sw):
        """Функция Баклея-Леверетта"""
        krw = self.relative_permeability_water(Sw)
        kro = self.relative_permeability_oil(Sw)

        # Добавляем малое число для избежания деления на ноль
        M = (krw / self.mu_water) / (kro / self.mu_oil + 1e-10)
        return M / (1 + M)

    def capillary_pressure(self, Sw):
        """
        Функция капиллярного давления по модели Брукса-Кори с плавным переходом
        в граничных зонах для повышения численной стабильности.
        """
        # Параметры Брукса-Кори (добавить как атрибуты класса в __init__)
        # self.entry_pressure = 0.3  # давление входа, МПа
        # self.pore_distribution_index = 2.0  # индекс распределения пор
        # self.wettability_factor = 0.8  # коэффициент смачиваемости (0-1)

        # Граничные значения насыщенности
        Swc = self.initial_water_saturation
        Sor = self.residual_oil_saturation

        # Избегаем численных проблем у границ диапазона насыщенности
        epsilon = 0.01  # параметр сглаживания вблизи границ

        if Sw <= Swc + epsilon:
            # Плавный переход к максимальному капиллярному давлению
            alpha = (Sw - Swc) / epsilon
            max_pc = self.entry_pressure * 3.0  # максимальное капиллярное давление, МПа
            return max_pc * (1.0 - alpha) + self.entry_pressure * alpha

        elif Sw >= 1 - Sor - epsilon:
            # Плавный переход к нулю капиллярного давления
            alpha = (1 - Sor - Sw) / epsilon
            return self.entry_pressure * 0.05 * alpha  # близко к нулю в конечной точке

        else:
            # Нормализованная водонасыщенность (эффективная)
            Se = (Sw - Swc) / (1 - Swc - Sor)

            # Модель Брукса-Кори
            pc = self.entry_pressure * (Se ** (-1.0 / self.pore_distribution_index))

            # Корректировка с учетом смачиваемости
            # Для гидрофобной среды (oil-wet) увеличиваем капиллярное давление
            pc = pc * (2.0 - self.wettability_factor)

            return pc

    def diffusion_coefficient(self, Sw):
        """Коэффициент капиллярной диффузии"""
        # Предотвращаем выход за граничные значения
        Sw = max(min(Sw, 0.99), 0.01)

        # Вычисление производной функции Баклея-Леверетта
        delta = 1e-4
        Sw_minus = max(Sw - delta, 0.01)
        Sw_plus = min(Sw + delta, 0.99)

        df_dS = (self.fractional_flow(Sw_plus) - self.fractional_flow(Sw_minus)) / (2 * delta)

        # Вычисление производной капиллярного давления
        dpc_dS = (self.capillary_pressure(Sw_plus) - self.capillary_pressure(Sw_minus)) / (2 * delta)

        # Увеличиваем коэффициент проницаемости для усиления эффекта
        k = 1.0  # изменить с 0.1 на 1.0

        mu = max(self.mu_water * Sw + self.mu_oil * (1 - Sw), 0.1)

        # Теоретическая формула с БОЛЬШИМ усилением
        D = -k / (self.porosity * mu) * df_dS * dpc_dS

        # Увеличиваем множитель с 0.05 до 1.0 для усиления эффекта
        max_diffusion = 0.45 * self.dx ** 2 / self.dt

        # Обрабатываем отрицательные значения и применяем ограничение устойчивости
        if D < 0:
            return min(abs(D) * 1.0, max_diffusion)  # изменить множитель с 0.05 на 1.0
        else:
            return min(D * 1.0, max_diffusion)  # изменить множитель с 0.05 на 1.0

    def run_simulation(self):
        """Запуск моделирования"""
        # Моделирование с учетом капиллярных эффектов
        for n in range(self.nt - 1):
            for i in range(1, self.nx):
                # Апвинд схема для конвективного члена
                f_i = self.fractional_flow(self.Sw_with_cap[n, i])
                f_im1 = self.fractional_flow(self.Sw_with_cap[n, i - 1])

                # Диффузионный член (капиллярные эффекты)
                D_i = self.diffusion_coefficient(self.Sw_with_cap[n, i])

                # Схема апвинд с учетом капиллярных эффектов
                self.Sw_with_cap[n + 1, i] = self.Sw_with_cap[n, i] - \
                                             (self.dt / self.dx) * (f_i - f_im1) + \
                                             (self.dt / self.dx ** 2) * D_i * (
                                                         self.Sw_with_cap[n, i + 1] - 2 * self.Sw_with_cap[n, i] +
                                                         self.Sw_with_cap[n, i - 1])

            # Граничное условие на правом конце
            self.Sw_with_cap[n + 1, -1] = self.Sw_with_cap[n + 1, -2]

        # Моделирование без учета капиллярных эффектов
        for n in range(self.nt - 1):
            for i in range(1, self.nx):
                # Апвинд схема для конвективного члена
                f_i = self.fractional_flow(self.Sw_without_cap[n, i])
                f_im1 = self.fractional_flow(self.Sw_without_cap[n, i - 1])

                # Схема апвинд без учета капиллярных эффектов
                self.Sw_without_cap[n + 1, i] = self.Sw_without_cap[n, i] - \
                                                (self.dt / self.dx) * (f_i - f_im1)

            # Граничное условие на правом конце
            self.Sw_without_cap[n + 1, -1] = self.Sw_without_cap[n + 1, -2]

    def calculate_recovery_factor(self):
        """Расчет коэффициента нефтеотдачи"""
        initial_oil = 1 - self.initial_water_saturation

        recovery_with_cap = np.zeros(self.nt)
        recovery_without_cap = np.zeros(self.nt)

        for n in range(self.nt):
            # Средняя нефтенасыщенность
            avg_oil_with_cap = 1 - np.mean(self.Sw_with_cap[n, :])
            avg_oil_without_cap = 1 - np.mean(self.Sw_without_cap[n, :])

            # Коэффициент нефтеотдачи
            recovery_with_cap[n] = (initial_oil - avg_oil_with_cap) / initial_oil
            recovery_without_cap[n] = (initial_oil - avg_oil_without_cap) / initial_oil

        return recovery_with_cap, recovery_without_cap

    def get_breakthrough_time(self):
        """Определение времени прорыва воды"""
        threshold = self.initial_water_saturation + 0.05

        # Время прорыва с учетом капиллярных эффектов
        breakthrough_with_cap = self.days
        for n in range(self.nt):
            if self.Sw_with_cap[n, -1] > threshold:
                breakthrough_with_cap = self.t[n]
                break

        # Время прорыва без учета капиллярных эффектов
        breakthrough_without_cap = self.days
        for n in range(self.nt):
            if self.Sw_without_cap[n, -1] > threshold:
                breakthrough_without_cap = self.t[n]
                break

        return breakthrough_with_cap, breakthrough_without_cap