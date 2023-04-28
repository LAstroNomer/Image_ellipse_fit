# Image_ellipse_fit
 Программа для определения неэллиптических объектов на астрономичесских исзображениях
 
 Данная программа реализуется в следующей логической последовательности:
 1. На вход подаётся изображение с вычтенным фоном. По нему,используя алгоритм 3sigmaclipping'а создаётся маска объектов.
 2. Данные объекты разделяются на отдельные источники.
 3. По молученным маскам создаются контуры и выгружаются в отдельный json файл
 4. Каждый эллипс приюлижается эллипсом (см. подробности ниже)
 5. По невязкам оценивается близость контура к эллиптическому. Создаётся текстовый документ с определёнными параметрами.
 
 ## Необходимые пакеты

 Для установки необходимых пакетов используйте 
 
 $pip3 install json astropy photutils 
 
 ## Тестовый запуск
 
 Для провеки работы программы запустите файл example.py
 
 $python3 example.py
 
 ## Запуск программы

 Для запуска программы используйте
 
 run.py --input INPUT.fits
 
 Для более тонкой настройки программы см. ниже пункт "файл настройки"
 
