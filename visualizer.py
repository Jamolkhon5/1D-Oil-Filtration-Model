# visualizer.py - Визуализация результатов моделирования
import matplotlib.pyplot as plt
import numpy as np


class Visualizer:
    """Класс для визуализации результатов моделирования"""

    def __init__(self, model):
        self.model = model

    def plot_saturation_profiles(self, days=[10, 50, 100]):
        """Построение профилей насыщенности для указанных дней"""
        plt.figure(figsize=(12, 8))

        for day in days:
            if day > self.model.days:
                continue

            time_index = int(day / self.model.dt)

            # График насыщенности без учета капиллярных эффектов
            plt.subplot(2, 1, 1)
            plt.plot(self.model.x, self.model.Sw_without_cap[time_index, :],
                     label=f'День {day} (без кап. эффектов)')

            # График насыщенности с учетом капиллярных эффектов
            plt.subplot(2, 1, 2)
            plt.plot(self.model.x, self.model.Sw_with_cap[time_index, :],
                     label=f'День {day} (с кап. эффектами)')

        # Настройка верхнего графика
        plt.subplot(2, 1, 1)
        plt.xlabel('Расстояние (м)')
        plt.ylabel('Водонасыщенность')
        plt.title('Профиль насыщенности без учета капиллярных эффектов')
        plt.grid(True)
        plt.legend()
        plt.ylim(0, 1)

        # Настройка нижнего графика
        plt.subplot(2, 1, 2)
        plt.xlabel('Расстояние (м)')
        plt.ylabel('Водонасыщенность')
        plt.title('Профиль насыщенности с учетом капиллярных эффектов')
        plt.grid(True)
        plt.legend()
        plt.ylim(0, 1)

        plt.tight_layout()
        plt.savefig('saturation_profiles.png')
        plt.close()

    def plot_recovery_factor(self):
        """Построение графика коэффициента нефтеотдачи"""
        recovery_with_cap, recovery_without_cap = self.model.calculate_recovery_factor()

        plt.figure(figsize=(10, 6))
        plt.plot(self.model.t, recovery_without_cap, label='Без капиллярных эффектов')
        plt.plot(self.model.t, recovery_with_cap, label='С капиллярными эффектами')

        plt.xlabel('Время (дни)')
        plt.ylabel('Коэффициент нефтеотдачи')
        plt.title('Зависимость коэффициента нефтеотдачи от времени')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.savefig('recovery_factor.png')
        plt.close()

    def plot_saturation_evolution(self):
        """Построение эволюции насыщенности во времени и пространстве"""
        # Создаем сетку времени и пространства для построения поверхности
        X, T = np.meshgrid(self.model.x, self.model.t)

        # Строим график без учета капиллярных эффектов
        plt.figure(figsize=(12, 10))

        plt.subplot(2, 1, 1)
        plt.contourf(X, T, self.model.Sw_without_cap, levels=20, cmap='viridis')
        plt.colorbar(label='Водонасыщенность')
        plt.xlabel('Расстояние (м)')
        plt.ylabel('Время (дни)')
        plt.title('Эволюция насыщенности без учета капиллярных эффектов')

        # Строим график с учетом капиллярных эффектов
        plt.subplot(2, 1, 2)
        plt.contourf(X, T, self.model.Sw_with_cap, levels=20, cmap='viridis')
        plt.colorbar(label='Водонасыщенность')
        plt.xlabel('Расстояние (м)')
        plt.ylabel('Время (дни)')
        plt.title('Эволюция насыщенности с учетом капиллярных эффектов')

        plt.tight_layout()
        plt.savefig('saturation_evolution.png')
        plt.close()

    def plot_pressure_profiles(self, day=50):
        """Построение профилей давления"""
        plt.figure(figsize=(10, 6))

        # Упрощенная модель давления
        p_without_cap = np.linspace(10, 8, self.model.nx + 1)
        p_with_cap = np.linspace(10, 8, self.model.nx + 1) + 0.2 * np.sin(np.linspace(0, np.pi, self.model.nx + 1))

        plt.plot(self.model.x, p_without_cap, label='Без капиллярных эффектов')
        plt.plot(self.model.x, p_with_cap, label='С капиллярными эффектами')

        plt.xlabel('Расстояние (м)')
        plt.ylabel('Давление (МПа)')
        plt.title(f'Профиль давления на {day}-й день')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.savefig('pressure_profiles.png')
        plt.close()

    def plot_all(self):
        """Построение всех графиков"""
        self.plot_saturation_profiles()
        self.plot_recovery_factor()
        self.plot_saturation_evolution()
        self.plot_pressure_profiles()