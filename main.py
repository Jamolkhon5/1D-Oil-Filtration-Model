# main.py - Главный файл для запуска моделирования
from model import OilFiltrationModel
from console_output import ConsoleOutput
from visualizer import Visualizer


def main():
    print("Запуск моделирования одномерной фильтрации нефти в пористой среде")
    print("с использованием метода апвинд")
    print("-" * 70)

    # Создаем модель
    model = OilFiltrationModel()

    print("Параметры модели:")
    print(f"Длина пласта: {model.length} м")
    print(f"Пористость: {model.porosity}")
    print(f"Вязкость нефти: {model.mu_oil} мПа·с")
    print(f"Вязкость воды: {model.mu_water} мПа·с")
    print(f"Начальная водонасыщенность: {model.initial_water_saturation}")
    print(f"Остаточная нефтенасыщенность: {model.residual_oil_saturation}")
    print(f"Количество узлов сетки: {model.nx}")
    print(f"Временной шаг: {model.dt} дней")
    print(f"Общее время моделирования: {model.days} дней")
    print("-" * 70)

    print("Запуск моделирования...")
    # Запускаем расчет
    model.run_simulation()
    print("Моделирование завершено!")
    print("-" * 70)

    # Создаем объекты для вывода и визуализации
    console = ConsoleOutput(model)
    visualizer = Visualizer(model)

    # Выводим результаты в консоль
    console.print_saturation_profile(day=50)
    console.print_recovery_factor()
    console.print_front_parameters()
    console.print_pressure_distribution(day=50)

    # Визуализируем результаты
    print("\nСоздание графиков...")
    visualizer.plot_all()
    print("Графики сохранены в текущей директории.")


if __name__ == "__main__":
    main()