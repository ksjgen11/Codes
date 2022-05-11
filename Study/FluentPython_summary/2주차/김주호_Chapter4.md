# Chap4. 텍스트와 바이트 (요약먼저)

  '1문자 == 1바이트'라는 개념을 거부하면서 4장을 시작한다. 전 세계적으로 유니코드가 채택되면서(웹사이트의 80%는 이미 UTF-8을 사용하고 있다) **텍스트 문자열**이라는 개념을 파일에 저장된 내용을 나타내는 **이진 시퀀스**와 분리해야 했으며, 파이썬 3는 문자열과 이진 시퀀스의 분리를 요구한다.

  **bytes, butearray, memoryview** 등의 이진 시퀀스 자료형을 간략히 살펴본 후 인코딩과 디코딩에 대해 살펴봤다. 그리고 몇 가지 코덱을 설명한 수 파이썬 소스 파일을 잘못 인코딩했을 때 발생하는 UnicodeEncodeError, UnicodeDecodeError, SyntaxError를 예방하거나 처리하는 방법을 설명한다.

  코드를 유지보수할 사람이 비아스키 문자로 구성된 자연어를 사용하고자 한다면 

- 파이썬2에서 코드를 실행할 필요가 없다면 

자연어를 식별자로 사용하는 것도 좋다. 그렇지만 전 세계 사람이 코드 개발에 참여한다면 아스키 문자로 구성된 영문단어로 식별자를 만들어야 한다.

  그러고 나서 메타데이터가 없을 때, 인코딩 방식을 탐지하는 이론과 방법을 살펴보았다. 이론적으로는 불가능하지만 실제로 **Chardet 패키지**는 여러 주류 인코딩 방식에 대해 텍스트의 코덱을 상당히 잘 찾아낸다. 그리고 UTF-16과 UTF-32에서 (그리고 UTF-8에서도 종종) 인코딩 방식을 알려주기 위해 사용하는 바이트 순서 표시를 설명한다.



  그리고 텍스트 파일 열기 연산에 대해 설명한다. 텍스트 파일을 열 때 **encoding 인수가 필수 인수는 아니지만, 반드시 사용하는 것이 좋다.** 인코딩 방식을 지정하지 않으면 프로그램은 기본 인코딩 방식이 여러 플랫폼에서 호환되지 않는 '평문'을 생성하게 된다. 그리고 파이썬이 기본 인코딩 방식을 알아내기 위해 사용하는 **locale.getpreferredendoding(), sys.getfilesystemencoding(), sys.getdefaultencoding()** 메서드 및 표준 입출력 파일에 대한 인코딩(sys.stdout.encoding emd)에 대해서 설명한다. 윈도우의 경우 안타깝게도 동일 컴퓨터 안에서도 이 설정들이 서로 다른 값을 가지며, 서로 호환되지 않는 경우도 있다. 반면 GNU/리눅스 및 OSX의 경우에는 거의 모든 곳에서 UTF-8을 사용한다.

  유니코드는 동일 문자를 여러 방식으로 표현할 수 있으므로 문자열 비교가 의외로 복잡하다. 따라서 텍스트를 매칭하려면 반드시 **문자열을 정규화**해야 한다. 정규화와 케이스 폴딩 외에도 악센트를 모두 제거하는 방법 등 필요에 따라 텍스트를 상당히 많이 변환하는 유틸리티 함수도 몇 가지 살펴본다. 그러고 나서 표준 **locale 모듈**을 활용해서 **유니코드 텍스트를 제대로 정렬하는 방법** 및 주의할 점을 알아보았다. 또한 까다로운 지역 설정에 의존하지 않고 외부 PyUCA 패키지를 사용해서 정렬하는 방법도 설명한다.

  마지막으로 유니코드 데이터베이스 및 이중 모드 API에 대해 간략히 설명했다. 유니코드 데이터베이스는 모든 문자에 대한 메타데이터를 제공하며, 이중모드 API는 re 및 os 모듈처럼 str이나 bytes 인수를 모두 받을 수 있으며, 인수의 종류에 따라 적절히 처리하는 함수를 제공한다.

## 4.2 바이트에 대한 기본 지식

이진 시퀀스를 위해 사용되는 내장 자료형은 

- bytes : 불변형
- bytearray : 가변형

2가지가 있다.



##  4.3 기본 인코더/디코더

텍스트를 바이트로 혹은 바이트를 텍스트로 변환하기 위해 파이썬 배포본에는 100여개의 코덱(인코더/디코더)이 포함되어 있다. 각 코덱은 utf_8과 같은 이름을 갖고 있는데, utf8, utf-8, U8 등으로 불리기도 한다. 코덱은 open(), str.encode(), butes.decode() 등의 함수를 호출할 때 encoding 인수에 전달해서 사용할 수 있다.


```python
for codec in ['latin_1', 'utf_8', 'utf_16']:
    print(codec, 'El Niño'.encode(codec), sep='\t')
```

    latin_1	b'El Ni\xf1o'
    utf_8	b'El Ni\xc3\xb1o'
    utf_16	b'\xff\xfeE\x00l\x00 \x00N\x00i\x00\xf1\x00o\x00'



(4.4 - 4.9는 코드예시를 정리하기 보다  keyword를 정리하고, search 하는 방향으로 학습.)

## 4.4 인코딩 / 디코딩 문제 이해하기

Although there is a generic UnicodeError exception, almost always the error reported is more specific: either an UnicodeEncodeError, when converting str to binary se‐quences or an UnicodeDecodeError when reading binary sequences into str. Loading Python modules may also generate SyntaxError when the source encoding is unex‐pected. We’ll show how to handle all of these errors in the next sections.

## 4.5 텍스트 파일 다루기
The best practice for handling text is the “Unicode sandwich” This means that bytes should be decoded to str as early as possible on input, e.g. when opening a file for reading. The “meat” of the sandwich is the business logic of your program, where text handling is done exclusively on str objects. You should never be encoding or de‐coding in the middle of other processing. On output, the str are encoded to bytes as late as possible. Most Web frameworks work like that, and we rarely touch bytes when using them. In Django, for example, your views should output Unicode str; Django itself takes care of encoding the response to bytes, using UTF-8 by default.

## 4.6 제대로 비교하기 위해 유니코드 정규화하기
String comparisons are complicated by the fact that Unicode has combining characters: diacritics and other marks that attach to the preceding character, appearing as one when printed.
For example, the word “café” may be composed in two ways, using 4 or 5 code points, but the result looks exactly the same.

## 4.7 유니코드 텍스트 정렬하기
Python sorts sequences of any type by comparing the items in each sequence one by one. For strings, this means comparing the code points. Unfortunately, this produces unacceptable results for anyone who uses non-ASCII characters.

## 4.8 유니코드 데이터베이스
The Unicode standard provides an entire database — in the form of numerous struc‐ tured text files — that includes not only the table mapping code points to character names, but also lot of metadata about the individual characters and how they are related. For example, the Unicode database records whether a character is printable, is a letter, is a decimal digit or is some other numeric symbol. That’s how the str methods isidentifier, isprintable, isdecimal and isnumeric work. str.casefold also uses in‐ formation from a Unicode table.


## 4.9 이중모드 str 및 bytes API
The standard library has functions that accept str or bytes arguments and behave differently depending on the type. Some examples are in the re and os modules.