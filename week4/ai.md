# Model Setup Details
The model used to help in completing this deliverable is OpenAI's gpt-5 (high reasoning mode) used in Agent mode on Cursor.

# Prompt 1 (Gemini 2.5 Pro)

## User Message
```
**Goal:** Based on the assignment description and Codon's context below, create an EPIC with Tickets to complete the assignment that asks to implement these algorithms in Codon and Python. 

**Assignment description:** 
```
Week 4 deliverable
The goal of this week’s deliverable is to learn how to implement various edit distance and dynamic programming algorithms.
Your job is to implement the following algorithms with the Codon:

global alignment;
local alignment;
semi-global (fitting) alignment; and
affine gap penalty global alignment.
You should implement both Codon and Python versions (which should be identical; it’s fine if they are not, though). All alignments should be tested with the following tests:

MT_human vs. MT_orang;
All pairs q1 vs. t1 … q5 vs. t5
The content of these sequences is to be found at https://github.com/lh3/ksw2/blob/master.
Use +3 for match, -3 for mismatch and -2 for gaps except for the affine alignment where you should use -5 for gap open and -1 for gap extension.

What (and when) to submit?
Submit a single link (the repository URL pointing to the commit to be graded) via Brightspace before Tue, Oct 21, at 23:59. The final commit must be successfully checked by the GitHub CI.
This assignment should be done alone. Use the same setup as in Week 1.
Your commit must have a working evaluate.sh (or any other extension) that is called from CI and that produces the following output:

Method            Language    Runtime
--------------------------------------global-mt_human   python      2000msglobal-q1         python      2000ms
local-q2          codon       1000msrun code snippet
Visit Manage Class to disable runnable code snippets×
The runtime is the time needed to run all tests (in milliseconds).

Grading
You will get one point for each method. Please use the following CI setup: https://github.com/inumanag/fall25-csc-bioinf/blob/main/.github/workflows/actions.yml.
You are expected to understand the code. You might be called to defend your submission at random during office hours.
Please describe all steps and gotchas in report.md in week4/ directory. Also include the time needed to complete this assignment (an estimate is OK). And don’t forget ai.md as well.
porting
0

30 views
6 Followup Discussions
Resolved
@65_f1

Anonymous Scale  6 days ago
Hi Professor, 
The due date of this deliverable is Sunday October 21st. However, the 21st is a Tuesday. Do you mean Sunday October 19? 
Thank you
0

Ibrahim Numanagić  6 days ago
Fixed!
1
Reply
Resolved
@65_f2

Wenhui Ma  2 days ago
Hi Professor,

I couldn’t find a definition of semi-global alignment in the lecture slides.
Just to confirm, does semi-global alignment require the entire shorter sequence to be aligned, or is it enough that one sequence reaches its end?

Thanks!
0

Ibrahim Numanagić  2 days ago
It is also known as fitting alignment. That is mentioned in the BWT slides.
0
Reply
Resolved
@65_f3

Ashley Osuna  2 days ago
"except for the semi-global alignment where you should use -5 for gap open and -1 for gap extension" should this be for affine gap penalty global alignment?
1

Ibrahim Numanagić  2 days ago
Yes, affine—thanks for correction!
1
```

**Codon context:** 

Codon is a Python compiler that type-checks Python code ahead-of-time and compiles it to native code.
Codon places some restrictions on Python, such as requiring all types to be knowable ahead of time, and requiring lists and other collections to have consistent types.
A non-homogeneous collection in Python must be converted to a legal equivalent in Codon. Here are some legal equivalents: firstly, the list can be converted to a list of tuples (example: Python type: List[List[int, str]], legal equivalent in Codon: List[Tuple[int, str]]).
Secondly, the list can be converted to a dictionary(example: Python type: List[List[int, str]]), legal equivalent in Codon: dict[int, str]), and so on.
The keywords: List, Tuple, and Dict must be used for annotating a list, a tuple and a dictionary, respectively and only when necessary.
Unless the element type can be inferred later by Codon, lists, sets and dictionaries cannot be initialized to empty literals without an annotation of the collection type (examples: list1 = [] by itself is not allowed. list1: List[int] = [] is allowed. list1 = [] is allowed if there is a list1.append(42) later in the code, because Codon will infer that the list element type is int).
There is no need to account for all of a list’s members types in its annotation. For instance a list that contains 3 ints is annotated as List[int]; using List[int, int, int] is illegal in Codon.
Also, Any is not a type in Codon and no variable can have the type Any. Using type Any is illegal.
Collection types must specify their element types in Codon: List or list alone are illegal, but List[int] or List[float] are legal.
In Codon, the types float and int cannot be used interchangeably and care should be taken when initializing variables to an integer if that variable should really be a float. Example: a = 10; a /= 3 is illegal in Codon and should be replaced by a = 10.0; a /= 3 to compile.
The same goes for collections of integers: a = [0 for _ in range(10)]; a[i] = 1.5 is illegal in Codon and should be replaced by a = [0. for _ in range(10)]; a[i] = 1.5. List[int] and List[float] are not compatible in Codon: a: List[float] = [1, 2, 3] is illegal in Codon because the right-hand side is first parsed as a list of integers. b: List[int] = [1.2, 3.4] is illegal in Codon because a list of floats cannot be assigned to a variable of type List[int].
Constant float values like inf, -inf, pi, etc. from the math library or any other library cannot be assigned to a type int unless they have been defined as a type int. Examples of statements that are illegal in Codon: m: int = math.inf, m: int = inf, m: int = float(‘inf’).
Classes in Codon require all members to be listed as fields along with their types, similar to Python’s dataclasses. Example: class A: def init(self, a): self.a = a is illegal in Codon and must be changed to class A: a: int def init(self, a): self.a = a.
All classes in Codon must list their fields in this way. A class that does not list its fields explicitly is illegal and will not compile in Codon. Constructor arguments themselves do not need to be explicitly typed in Codon.
Codon does not support assigning methods in a class definition: class A: def foo(self): return 42 bar = foo is illegal in Codon and must be changed to class A: def foo(self): pass def bar(self): return foo()
Variables or fields that can be assigned to None must be marked as optional. For example, this code is illegal: class Node: def init(self, left = None, right = None): self.left = left self.right = right and must be changed to class Node: left: Optional[Node] right: Optional[Node] def init(self, left: Optional[Node] = None, right: Optional[Node] = None): self.left = left self.right = right Codon does not support type names in quotes; such names must be given without quotes: foo: ‘Foo’ is illegal in Codon, and foo: Foo the valid replacement. Similarly: class Node: next: Optional[‘Node’] is illegal in Codon and must be written as class Node: next: Optional[Node]
Codon does not require type annotations on variables, function arguments or function return types, as it will infer these types automatically.
Do not add new type annotations on variables, function arguments or function return types in cases where these types can be inferred. Codon knows the types of all literals (e.g. 42, 3.14, etc.), meaning variables assigned to literals (as in x = 42) do not need to be type annotated.
Codon knows the return types of all the Python standard library functions (e.g. int(), time() etc.), meaning variables assigned to the result of these functions (as in x = int(3.14)) do not need to be type annotated.
Codon does not support string formatting via % or .format(): ‘%f’ % x is illegal in Codon and should be replaced with a format-string version like f’{x}'.
Codon includes a new NumPy implementation, so NumPy arrays and APIs can be used regularly in Codon. However, Codon’s ndarray type is parameterized by its data type and dimension: ndarray[float, 2] is a 2D array of floats, ndarray[int, 1] is a 1D array of int.
Codon deduces output type and dimension of most NumPy operations, unless it is impossible to do so as is the case with np.squeeze().
If the imported libraries are not implemented by Codon, you must modify the import. For instance ‘import requests’ must be transformed to ‘from python import requests’, but if the original import already starts with ‘from’ like ‘from bs4 import BeautifulSoup’, you have to change it to ‘from python import bs4’, and then use bs4.BeautifulSoup wherever required.   These were just some demonstrative samples to teach you the import rules.
math, random, time and other common libraries are all included in Codon and you do not need to import them from Python.
```

## LLM Response

```
Here is an EPIC with corresponding tickets to complete the assignment, based on the provided description and Codon's context.

EPIC: Implement and Benchmark Alignment Algorithms in Python and Codon
Description: Implement, test, and benchmark four dynamic programming alignment algorithms (Global, Local, Semi-Global, Affine Gap Global) in both standard Python and Codon. The project must be hosted in a GitHub repository with a CI pipeline that uses evaluate.sh to run all specified tests (MT_human vs. MT_orang, q1-t5 vs. t1-t5) and report performance in the required format. The final submission includes documentation on the process and AI usage.

Tickets
Ticket 1: Project Setup and Data Acquisition 📚
User Story: As a student, I need to set up my repository structure and acquire all necessary test data so I can begin development and testing.

Acceptance Criteria:

Create the week4/ directory within the course GitHub repository.

Create a data/ subdirectory inside week4/.

Download and place the following files into data/ using the content provided in the prompt:


MT-human.fa 


MT-orang.fa 


q1.fa (containing q1-q5) 


t1.fa (containing t1-t5) 

Create placeholder files: week4/report.md and week4/ai.md.

Initialize a Python virtual environment and a requirements.txt file (e.g., for numpy).

Ticket 2: Create a Multi-FASTA Parser 🧬
User Story: As a developer, I need a utility function to read and parse the FASTA files so I can easily extract sequences for the alignment algorithms.

Acceptance Criteria:

Create a function (e.g., parse_fasta(filename)) in a Python module.

The function must correctly read a FASTA file and return a dictionary or a list of tuples, mapping headers to their corresponding sequences (e.g., {'q1': 'AG', 'q2': 'AGCTAcCGCAT', ...}).

The parser must correctly handle multi-line sequences (like in MT-human.fa).

The parser must be compatible with both Python and Codon. (This may require type annotations like Dict[str, str] for Codon compliance).

Ticket 3: Python Implementation - Global Alignment (Needleman-Wunsch)
User Story: As a developer, I need a Python function that performs global alignment using the Needleman-Wunsch algorithm.

Acceptance Criteria:

Implement global_align(seq1: str, seq2: str) in a Python module (e.g., alignment.py).

The function must use a dynamic programming matrix (e.g., numpy.ndarray).

Uses the specified scoring: +3 for match, -3 for mismatch, and -2 for gaps.

The function must return the final alignment score.

Ticket 4: Python Implementation - Local Alignment (Smith-Waterman)
User Story: As a developer, I need a Python function that performs local alignment using the Smith-Waterman algorithm.

Acceptance Criteria:

Implement local_align(seq1: str, seq2: str) in the Python module.

The function must use a DP matrix where no cell value can be negative (i.e., max(..., 0)).

Uses the specified scoring: +3 for match, -3 for mismatch, and -2 for gaps.

The function must return the maximum score found anywhere in the matrix.

Ticket 5: Python Implementation - Semi-Global (Fitting) Alignment
User Story: As a developer, I need a Python function that performs semi-global (fitting) alignment.

Acceptance Criteria:

Implement semi_global_align(seq1: str, seq2: str) in the Python module.

The DP matrix initialization must not penalize starting gaps (e.g., the first row and first column are initialized to 0).

Uses the specified scoring: +3 for match, -3 for mismatch, and -2 for gaps.

The function must return the maximum score found in the last row or last column.

Ticket 6: Python Implementation - Affine Gap Global Alignment
User Story: As a developer, I need a Python function that performs global alignment using an affine gap penalty.

Acceptance Criteria:

Implement affine_global_align(seq1: str, seq2: str) in the Python module.

The implementation must use three DP matrices (e.g., M, I_x, I_y).

Uses the specified scoring: +3 for match, -3 for mismatch.

Uses the specified gap penalties: -5 for gap open and -1 for gap extension.

Returns the final optimal alignment score.

Ticket 7: Port Alignment Algorithms to Codon ⚡
User Story: As a developer, I need to port my Python alignment module to Codon to create a high-performance, compiled version.

Acceptance Criteria:

Create a Codon version of the alignment module (e.g., alignment.codon).

Modify the Python code to be Codon-compliant, paying attention to the Codon context.

Ensure all DP matrices (e.g., numpy.ndarray) are correctly typed (e.g., ndarray[int, 2]).

Ensure any use of infinity (e.g., for affine gap initialization) uses float values (e.g., -math.inf or -float('inf')) and not assigned to int types.

Add type annotations for function arguments, return types, and collections only where Codon cannot infer them (e.g., List[Tuple[str, str]] for the FASTA parser result).

The alignment.codon file successfully compiles using codon build -release.

Ticket 8: Create Evaluation and Benchmarking Script (evaluate.sh) ⏱️
User Story: As a student, I need a shell script that runs all required algorithm/language/test-pair combinations and outputs the runtime in the specified format for grading.

Acceptance Criteria:

Create an executable week4/evaluate.sh script.

The script must compile the Codon code (if not already compiled).

The script must run all 4 algorithms in both Python and Codon versions.

Each algorithm/language pair must be tested on:

MT_human vs. MT_orang

q1 vs. t1

q2 vs. t2

q3 vs. t3

q4 vs. t4

q5 vs. t5

The script must time each individual run in milliseconds.

The script's standard output must exactly match the specified format (including the header):

Method          Language    Runtime
--------------------------------------
global-mt_human python      [time]ms
global-q1       python      [time]ms
...
local-q2        codon       [time]ms
...
Ticket 9: Configure GitHub Actions CI 🤖
User Story: As a student, I need my GitHub repository to automatically run the evaluation script via CI to verify my submission is working correctly.

Acceptance Criteria:

Create a .github/workflows/actions.yml file.

Use the CI setup file provided in the assignment: https://github.com/inumanag/fall25-csc-bioinf/blob/main/.github/workflows/actions.yml.

The CI workflow must successfully install Codon, install dependencies, and execute the week4/evaluate.sh script.

The final commit to be graded must show a "green check" (passing CI) in GitHub.

Ticket 10: Write Final Report and Documentation ✍️
User Story: As a student, I need to document my implementation, report any challenges, and declare my AI usage to complete the assignment.

Acceptance Criteria:

Complete week4/report.md.

The report must describe the steps taken and any "gotchas" encountered, especially during the Codon porting process.

The report must include an estimate of the time needed to complete the assignment.

Complete week4/ai.md to describe all uses of AI in completing the assignment.
```

# Prompt 2 (gpt-5-high agent mode)

## User Message

```
Ticket 2: Create a Multi-FASTA Parser 🧬
User Story: As a developer, I need a utility function to read and parse the FASTA files so I can easily extract sequences for the alignment algorithms.

Acceptance Criteria:

Create a function (e.g., parse_fasta(filename)) in a Python module.

The function must correctly read a FASTA file and return a dictionary or a list of tuples, mapping headers to their corresponding sequences (e.g., {'q1': 'AG', 'q2': 'AGCTAcCGCAT', ...}).

The parser must correctly handle multi-line sequences (like in MT-human.fa).

The parser must be compatible with both Python and Codon. (This may require type annotations like Dict[str, str] for Codon compliance).
```

## Agent Response

```
Thought for 40s

I'll first confirm Codon-compatible typing for returning a dictionary of sequences, then implement parse_fasta under week4/python/ to handle multi-line FASTA and mixed-case sequences.

--------------------------------------------------------------------

Called Codon's MCP to `get_similar_docs()` with the following params: 

{
  "query": "Does Codon support Python type annotations using typing.Dict vs built-in generics like dict[str, str]? What is recommended for Codon compatibility?",
  "k": 5
}

Result: 
{"result":["While Codon supports nearly all of Python's syntax, it is not a drop-in replacement, and large codebases might require modifications to be run through the Codon compiler. For example, some of Python's modules are not yet implemented within Codon, and a few of Python's dynamic features are disallowed. The Codon compiler produces detailed error messages to help identify and resolve any incompatibilities. Codon supports seamless Python interoperability to handle cases where specific Python libraries or dynamism are required, and also supports writing Python extension modules that can be imported and used from larger Python codebases.","Codon's type system is designed to be non-intrusive, meaning it will infer types as much as possible without manual type annotations. This way, Python programs typically \"just work\" without requiring code changes.","CPython? Codon tries to follow CPython's syntax, semantics and APIs as closely as possible, aside from a few cases where Codon differs from CPython for performance reasons (one example being Codon's 64-bit int vs. CPython's arbitrary- width int).","There are ways to use Codon in larger Python codebases via its JIT decorator or Python extension backend.","Codon supports Python classes as you would expect."]}

--------------------------------------------------------------------

Thought for 1m 21s

I'll create a new parser utility with a Codon-compatible signature, then run quick assertions against provided FASTA files to verify multi-line handling and header parsing.

Thought for 59s 

I'll implement the FASTA parser in a new week4/python/utils.py with Codon-compatible typing and then validate it against the provided FASTA files.

Created utils.py (+44 additions)

Thought for 28s

I'll refine typing in the parser to use Optional and List for Codon compatibility, then lint the file and run a quick Python one-liner to validate parsing of all provided FASTA files.

Made 3 changes in utils.py

I'll verify the implemented parser file content and run lint checks, then mark the implementation task complete and begin validation.

No linter errors found in utils.py

Thought for 30s 

I'll run the parser via the demo CLI in main.py against the four FASTA files to validate record counts and ensure multi-line concatenation works, then mark validation complete.

<c427/fall25-csc-bioinf/week4/python/data/MT-orang.fa
/Users/tarekalakkadp/Desktop/uvic/fourth-year/fall/csc427/fall25-csc-bioinf/week4/python/data/q1.fa -> 5 records
  q1: 2 bp
  q2: 11 bp
  q3: 54 bp
  q4: 4 bp
  q5: 34 bp
/Users/tarekalakkadp/Desktop/uvic/fourth-year/fall/csc427/fall25-csc-bioinf/week4/python/data/t1.fa -> 5 records
  t1: 3 bp
  t2: 19 bp
  t3: 72 bp
  t4: 15 bp
  t5: 34 bp
/Users/tarekalakkadp/Desktop/uvic/fourth-year/fall/csc427/fall25-csc-bioinf/week4/python/data/MT-human.fa -> 1 records
  MT_human: 16569 bp
/Users/tarekalakkadp/Desktop/uvic/fourth-year/fall/csc427/fall25-csc-bioinf/week4/python/data/MT-orang.fa -> 1 records
  MT_orang: 16499 bp
tarekalakkadp@MacBookPro fall25-csc-bioinf % 


Added a Codon-compatible FASTA parser and validated it.
- Implemented parse_fasta in week4/python/utils.py with Dict[str, str] return type, proper handling of multi-line sequences, and identifier extraction from headers.
- Added a small CLI demo in week4/python/main.py to print record counts and lengths for given FASTA files.
- Verified on q1.fa, t1.fa, MT-human.fa, MT-orang.fa: correct record counts and concatenated multi-line sequences.
Run:
python3 week4/python/main.py week4/python/data/q1.fa week4/python/data/t1.fa week4/python/data/MT-human.fa week4/python/data/MT-orang.fa

Files touched: week4/python/utils.py, week4/python/main.py
```
