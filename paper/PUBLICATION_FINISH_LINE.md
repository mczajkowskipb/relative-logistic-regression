# CRIT: zamknięta droga do publikacji (PAPER MODE)

**Decyzja: 26 października 2026.** Od tego momentu priorytetem jest manuskrypt, a nie kolejne eksperymenty eksploracyjne H5–H10. Wyniki negatywne publikujemy uczciwie; nie dostrajamy parametrów pod obejrzane dane.

## Stan rzeczy na wejściu

- Istnieje oddzielny zamrożony manuskrypt podstawowego RLR: Relative Logistic Regression: A Probabilistic Refinement of Hard-Voting Pair-Based Omics Classification. Nowy artykuł musi mieć odmienny wkład.
- H0/G: po udostępnieniu wszystkich 496 par dla 32 genów odnaleziono 12/12 prawdziwych struktur XOR (5 orientacji bezpośrednich, 7 zamienionych); to dane syntetyczne wykorzystane przy rozwoju.
- H3-fix1: ślepe wyszukiwanie dwóch relacji XOR bez oracle — 19/20 trafień top6 przy silnym syntetycznym sygnale, dane tylko D.
- H4A: nowy, zamrożony test 288 przypadków D. Trafienia top6: siła 2.4: 70/72; 1.6: 36/72; 1.0: 4/72. Łącznie 1728 zweryfikowanych wyników kandydatów.
- Dotąd NIE wykazano niezależnej skuteczności CRIT+LSH na testowym E, poprawy względem oryginalnego RLR na nowych kohortach biomedycznych ani kontroli błędu fałszywej akceptacji.

## Bramka P1 — jeden końcowy protokół D/V/E

**Jeden z góry ustalony eksperyment, nie wielokrotne poprawianie aż pojawi się pozytywny wynik.**

1. Zachować bez zmian oryginalne RLR, ich konfiguracje, obliczenia strict OOF i schemat D/V/E z v9. Wykazać zgodność predykcji ORIGINAL_RLR z zamrożonymi referencjami.
2. LSH z H3 (18 bitów, 24 tablice) generuje maksymalnie sześć wspólnych propozycji relacji wyłącznie na D. Wyniki H4A nie służą do dostrajania.
3. Zdefiniować jawnie i zamrozić mapowanie par XOR na konstrukcję CONDITIONAL_TREE oraz GLOBAL_INTERACTION, z zamiennością ról obu relacji. Wybór tylko jednego modelu w D; żadnego uczenia nominacji na V lub E.
4. V jednorazowo akceptuje albo odrzuca korektę. Na nietkniętym E mierzyć logloss jako pierwotny endpoint, Brier, predykcje na poziomie pacjenta, odsetek zaakceptowanych korekt i dokładną zgodność fallback z oryginalnym RLR.
5. Uwzględnić niezależne nowe syntetyczne scenariusze słaby/umiarkowany/silny XOR oraz null i uczciwie zliczyć wszystkie porażki. Ponieważ sam test H4A wykazał słabą wykrywalność, nie zakładać sukcesu.
6. W ramach tego samego zamrożonego pipeline przeprowadzić możliwie duży zestaw niezależnych danych null, docelowo kilkaset, jeśli budżet obliczeniowy pozwala. Przed twierdzeniem o kontroli alfa obliczyć właściwe przedziały niepewności i uwzględnić adaptacyjną selekcję oraz liczbę porównań. Jeśli próba jest zbyt mała, zrezygnować z tezy o kontroli błędu.
7. Niezależny audyt: brak przecieków pacjentów/etykiet, zgoda źródeł, pełna lista modeli, SHA, replikowalny split, predykcja każdego pacjenta i test fallback.
8. Bramka interpretacyjna: jeśli nie ma uczciwej poprawy na E, STOP — nie poprawiamy kolejny raz algorytmu na tych samych danych. Wybieramy artykuł o mechanizmie odkrywania i jego ograniczeniach.

## Bramka P2 — jeden niezależny benchmark biomedyczny, tylko jeśli P1 uzasadnia

- Wybrać przed rozpoczęciem **co najmniej dwie** dotychczas niewykorzystywane do tworzenia CRIT publiczne kohorty lub ściśle niewidziane grupy pacjentów.
- Zamrozić przetwarzanie, identyfikatory kohort, jednostkę walidacji, rozdzielenie pacjentów/grup, metrykę pierwotną E logloss, oryginalny RLR jako baseline i ramiona kontrolne.
- Raportować także wyniki negatywne, Brier, klasyfikację drugorzędną, stabilność wybranych par i koszt obliczeń. Nie traktować wielokrotnie obejrzanych Chen, Gordon, CompCancer i PTB jako nowej niezależnej walidacji.
- Jeśli nie ma danych lub zysku, nie formułować mocnej tezy o przewadze biomedycznej.

## Bramka P3 — redakcja publikacji zaczyna się TERAZ

- [x] Draft angielski z dowodami i ograniczeniami: paper/CRIT_manuscript_v0.md.
- [ ] Kontrola literatury: wcześniejsze drzewa oparte na parach (w tym Czajkowski/Krętowski 2011), logika kombinatoryczna/XOR, binary-signature LSH; sformułować wkład jako joint interaction discovery, nie wynalezienie drzew ani XOR.
- [ ] Jedna figura algorytmu (podzielić część wykonaną H3 od przyszłej integracji P1).
- [ ] Tabele z wynikami wszystkich 36 komórek H4A, oszacowaniami niepewności i SHA.
- [ ] Wyniki P1/P2 uzupełnić wyłącznie po rzeczywistym PASS i niezależnym sprawdzeniu.
- [ ] Ostateczny wybór czasopisma, komplet bibliografii, autorstwo i poprawna wersja LaTeX.
- [ ] Zamrozić źródła, numery, podpisy rycin i manuskrypt.

## Zakres odroczony, domyślnie nie wykonywać

- Dalsze H5/H6/H7 itd. jako otwarty program prób i korekt algorytmu.
- Zwiększanie liczby tablic/budżetów H3 według obejrzanych H4A, a następnie nazywanie tego niezależnym testem.
- GPU, wdrożenie na pełne transkryptomy, pakiety i szeroki benchmark przed weryfikacją zysku na E.

## Definicja końca

**Po P1 decydujemy o klasie artykułu:** jeżeli pełny model uzyskuje poprawę na uczciwym E, przechodzimy do P2 i celujemy w silniejszy artykuł empiryczny. W przeciwnym razie zapisujemy uczciwą publikację metodyczno-diagnostyczną o odkrywaniu interakcji (jeśli oryginalność i recenzje to uzasadnią), bez roszczeń do klinicznej przydatności.

Nie obiecujemy kalendarzowego terminu ukończenia: zależy on od rzeczywistych wyników P1, dostępności nowych kohort i oceny redakcyjnej. **Liczbę kolejnych otwartych etapów eksperymentalnych ograniczamy do dwóch.**
