## Compilação

```bash
cmake .
cmake --build .
```

## Utilização

```bash
./QOBI
```

Abre uma janela em que se pode desenhar compostos. Clicar nas telas C, O, N, S, F, L, B e I para alterar a espécie do átomo selecionado.

Carbono é representado a branco, oxigénio a vermelho, azoto a azul, enxofre a amarelo e os halogenos a verde. Ligações entre carbono e uma outra espécie ficam com a cor desta.

Para criar ligações duplas ou triplas, conectar outra vez os dois átomos já desenhados. Não há limite no número de ligações possíveis (os nomes gerados nestas situações não têm sentido).

Clicar no espaço alterna o modo de edição entre adição de átomos e o modo de movimento de átomos.

Clicar no enter para imprimir o nome da molécula no terminal. A cadeia principal fica destacada a rosa.

<img width="853" height="917" alt="image" src="https://github.com/user-attachments/assets/754a2128-b522-4f80-adc0-32f0a6c6607a" />

## Alcance

O algoritmo gera nomes de moléculas, segundo as regras de nomenclatura substitutiva da IUPAC.

Faz reconhecimento de vários grupos funcionais com oxigénio, azoto e enxofre (ácidos carboxílicos, ésteres, amidas, nitrilos, aldeídos, cetonas, álcoois, tióis, aminas e iminas). Não reconhece, por exemplo, ácidos sulfónios.

Também é possível desenhar halogenos (flúor, cloro, bromo e iodo). Não há suporte para outras espécies atómicas, nem moléculas carregadas. Ainda não há um sistema para nomear cadeias muito longas.

Há suporte de compostos cíclicos simples (incluindo benzeno), mas, neste momento, uma molécula policíclica origina uma falha no programa.
