find -name *.o | xargs rm -f
find -name *.gcno | xargs rm -f
find -name *.gcda | xargs rm -f
find -name *.gcov | xargs rm -f
find -name *.txt | xargs rm -f
rm libSYSV.* -f