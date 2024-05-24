# Метод Нелдера-Мида

## Введение

Метод Нелдера – Мида применяется для нахождения решения задачи оптимизации вещественных функций многих переменных

$$f(x) -> min, \forall x \in R^n$$

причем функция f(x) не является гладкой.

Другой особенностью метода является то, что на каждой итерации вычисляется значение функции f(x) не более чем в трех точках. Это особенно важно в случае сложно-вычислимой функции f(x). Метод Нелдера - Мида прост в реализации и полезен на практике, так как не подразумевает поиск решения через производные. Но, с другой стороны, для него <u>не существует теории сходимости*</u> - алгоритм может расходиться даже на гладких функциях.

Симплексом (n-симплексом) называется выпуклая оболочка множества** независимых (n + 1) точек (вершин симплекса).

> *Функция f(x) сходится к предельному значению L при x→a, если для любого положительного числа E существует положительное число G, такое что для всех x удовлетворяющих условиям 0\<\|x−a\|\<G, выполнено неравенство \|f(x)−L\|\<E.

> **Выпуклая оболочка множества - это наименьшее выпуклое множество, которое его содержит. Им будет пересечение всех выпуклых множеств, содержащих множество.

---

## Описание алгоритма

Алгоритм метода Нелдера – Мида включает следующие шаги:

 1. <u>Инициализация:</u> строится симплекс - многогранник в n-мерном пространстве, где n - количество переменных функции оптимизации;

 2. <u>Оценка значений функции:</u> вычисляются значения функции в вершинах симплекса;

 3. <u>Сортировка:</u> вершины сортируются в порядке убывания значений функции;

 4. <u>Центр тяжести:</u> вычисляется центр тяжести симплекса без наихудшей вершины;

 5. <u>Отражение:</u> проводится отражение на основе центра тяжести и лучшей вершины симплекса;

 6. <u>Оценка нового значения функции в отраженной точке;</u>

 7. <u>Расширение:</u> если новая точка лучше, чем лучшая вершина симплекса, происходит расширение в этом направлении;

 8. <u>Сжатие:</u> если новая точка хуже, чем наихудшая вершина, происходит сжатие симплекса в сторону лучшей вершины;

 9. <u>Метод средней точки:</u> если ни одно из вышеупомянутых действий не улучшает значение функции, то сжимаем симплекс в сторону центра тяжести без наихудшей вершины;

10. <u>Условие останова:</u> алгоритм завершается, если выполнено условие сходимости, например, при достижении заданной точности или после определенного числа итераций.

Алгоритм метода Нелдера – Мида повторяется до достижения оптимального значения функции или условия останова.

---

## Математические тесты

1. Функция Химмельблау — функция двух переменных, используемая для проверки эффективности алгоритмов оптимизации.

Она определяется формулой:

$$(x_1^2 + x_2 - 11)^2 + (x_1 + x_2^2 - 7)^2$$

Функция Химмельбау содержит четыре равнозначных локальных минимума:

$$f(3, 2)=0;$$

$$f(-2,80511…, 3,13131…)=0;$$

$$f(-3,77931…, -3,28318…)=0;$$

$$f(3,58442…, -1,84812…)=0.$$

#### Проведём несколько запусков алгоритма, изменяя коэффициенты

| Коэфф. отражения | Коэфф. растяжения | Коэфф. сжатия | Кол-во итераций | Значение x_1 | Значение x_2 | Значение функции |
|------------------|-------------------|---------------|-----------------|--------------|--------------|------------------|
| 0,5              | 2,0               | 0,5           | 143             | -2,81662876  | -1,88177226  | 63,86908948      |
| 1,0              | 2,0               | 0,5           | 80              | 2,99999595   | 1,99998769   | 0,00000000       |
| 1,5              | 2,0               | 0,5           | 143             | 5,33787090   | 4,58743657   | 863,21893117     |
| 1,0              | 1,5               | 0,5           | 72              | 2,99999822   | 2,00000891   | 0,00000000       |
| 1,0              | 2,0               | 0,5           | 80              | 2,99999595   | 1,99998769   | 0,00000000       |
| 1,0              | 2,5               | 0,5           | 96              | 2,99998766   | 2,00000149   | 0,00000001       |
| 1,0              | 2,0               | 0,3           | 56              | 3,04597535   | 1,88607480   | 0,18431711       |
| 1,0              | 2,0               | 0,5           | 80              | 2,99999595   | 1,99998769   | 0,00000000       |
| 1,0              | 2,0               | 0,7           | 120             | 3,00000769   | 2,00000796   | 0,00000000       |

Как видно из таблицы, при изменении коэффициентов изменяется количество итераций и точность ответа, при изменении коэффициента отражения ответ вообще может быть неправильным.

2. Функция Розенброка

Невыпуклая функция, используемая для оценки производительности алгоритмов оптимизации, предложенная Ховардом Розенброком в 1960 году. Считается, что поиск глобального минимума для данной функции является нетривиальной задачей. Является примером тестовой функции для локальных методов оптимизации. Имеет минимум 0 в точке (1,1).

Функция определена следующим образом:

$$(1 - x_1)^2 + 100*(x_2 - x_1^2)^2$$

#### Проведём несколько запусков алгоритма, изменяя коэффициенты

| Коэфф. отражения | Коэфф. растяжения | Коэфф. сжатия | Кол-во итераций | Значение x_1 | Значение x_2 | Значение функции |
|------------------|-------------------|---------------|-----------------|--------------|--------------|------------------|
| 0,5              | 2,0               | 0,5           | 144             | 0,07782825   | -0,00675289  | 0,86681067       |
| 1,0              | 2,0               | 0,5           | 156             | 1,00018018   | 1,00035833   | 0,00000003       |
| 1,5              | 2,0               | 0,5           | -               | -            | -            | -                |
| 1,0              | 1,5               | 0,5           | 272             | 0,99999559   | 0,99999178   | 0,00000000       |
| 1,0              | 2,0               | 0,5           | 156             | 1,00018018   | 1,00035833   | 0,00000003       |
| 1,0              | 2,5               | 0,5           | 195             | 1,00000299   | 1,00000239   | 0,00000000       |
| 1,0              | 2,0               | 0,3           | 107             | 0,99998876   | 0,99998032   | 0,00000000       |
| 1,0              | 2,0               | 0,5           | 156             | 1,00018018   | 1,00035833   | 0,00000003       |
| 1,0              | 2,0               | 0,7           | 509             | 1,00014612   | 1,00027995   | 0,00000004       |

Как видно из таблицы, при изменении коэффициентов изменяется количество итераций и точность ответа, при изменении коэффициента отражения ответ вообще может быть неправильным.

---

## Вывод

Метод Нелдера-Мида является эффективным инструментом для оптимизации функций без использования градиента. Его простота и эффективность делают его популярным выбором для различных задач оптимизации. Но данный метод не подходит для задач оптимизации со многими локальными минимумами.
