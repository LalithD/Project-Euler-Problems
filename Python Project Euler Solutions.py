## Problem 1: "Multiples of 3 and 5" Answer: 233168
def sumOfMultiples(multiple, max): # finds the sum of the multiples of parameter multiple from 0 to max (exclusive)
	import math
	numOfNumbers = math.ceil(max/multiple)-1
	return multiple * numOfNumbers * (numOfNumbers + 1) / 2 # uses formula 1+2+...+n = (n^2+n)/2
sumOfMultiples(3, 1000) + sumOfMultiples(5, 1000) - sumOfMultiples(15, 1000)

## Problem 2: "Even Fibonacci numbers" Answer: 4613732
def evenFib(n): # returns the sum of all even Fibonacci numbers less than n
	sum = 0
	a, b = 1, 0 # a is the -1th fibonacci number, b is the 0th fibonacci number
	while b < n: # after n iterations, b is the (3n)th fibonacci number, noticing that precisely (3n)th form the set of
		# even Fibonacci numbers
		a, b = b, a+b; a, b = b, a+b; a, b = b, a+b
		if b < n:
			sum += b
	return sum
evenFib(4 * 10**6)

## Problem 3: "Largest prime factor" Answer: 6857
def largestPrimeFactor(n): ## returns the largest prime factor of integer n
	import math
	if n % 2 == 0:
		return max(2, largestPrimeFactor(n/2))
	for i in range(3, math.ceil(n**0.5)+1, 2):
		if n % i == 0:
			return max(i, largestPrimeFactor(n/i))
	return n
largestPrimeFactor(600851475143)

## Problem 4: "Largest palindrome product" Answer: 90609
def largestPalindrome(n):
	maxPalindrome = 0
	for i in range(9 * 10**(n-1), 10**n): # assumes a palindrome exists in this range (self-evident)
		for j in range(9 * 10**(n-1), 10**n):
			if str(i * j) == str(i * j)[::-1]:
				maxPalindrome = max(maxPalindrome, i * j)
	return maxPalindrome
largestPalindrome(3)

## Problem 5: "Smallest multiple" Answer: 232792560
def smallestMultiple(n): # returns the smallest number that is a multiple of each number in array
	import fractions
	multiple = 1
	for i in range (1, n+1):
		multiple = multiple * i / fractions.gcd(multiple, i) # uses the property that lcm(a,b) * gcd(a,b) = a*b
	return multiple
smallestMultiple(20)
# this particular problem (for n=20) can also be solved given smallestMultiple(10) and multiplying by 11 * 13 * 17 * 19 * 2
# given a general smallestMultiple(n) and trying to find smallestMultiple(m), we want to multiply smallestMultiple(n) by all
# primes p for n < p <= m and all natural numbers q for which q^e <= n and n < q^(e+f) <= m, where f > 0. Note that We would multiply
# by q^f, for cases where f > 1
# This method is evident by noticing that lcm(a,b) * gcd(a,b) = a*b

## Problem 6: "Sum square difference" Answer: 25164150
def sumOfSquares(n): # returns the sum of the squares of the first n natural numbers
	return n * (n+1) * (2*n + 1) / 6
def sumOfIntegers(n): # returns the sum of the first n natural numbers
	return n * (n+1) / 2
def sumOfCumbes(n): # returns the sum of the cubes of the first n natural numbers (for completeness)
	return sumOfIntegers(n)**2
sumOfIntegers(100)**2 - sumOfSquares(100)
# Note: all of the above formulas can be proven by induction on n, noticing that 1^e+2^e+...+n^e can be expressed as
# an (e+1)th degree polynomial (again, we can prove this by induction on e)

## Problem 7: "10001st prime" Answer: 104743
def findNthPrime(n): # fast (28.5 seconds to find 100k prime)
	import math
	num = 3
	while n > 1:
		isPrime = True
		if num % 2 != 0:
			for i in range(3, math.floor(num**0.5)+2, 2):
				if num % i == 0:
					isPrime = False
					break
		else:
			isPrime = False
		if isPrime:
			n -= 1
		num += 2
	return num-2
findNthPrime(10001)

## Problem 8: "Largest Product in a series" Answer: 23514624000
def largestProduct(n, m): # will take a large integer n and find the largest product formed by m consecutive digits
	# goes through approximately n * m steps
	import math
	maximum = 1
	length = math.floor(math.log10(n))+1 # digits in n (length of n as a string)
	for i in range(0,length-m+1):
		s = str(n)[i:i+m] # contains an m element string of digits
		product = 1 # product will be the product of digits in s
		for j in range(0,m):
			product *= int(s[j])
		maximum = max(maximum, product) # replaces max if new string is larger
	return maximum
largestProduct(n, 13) # n is a specific 1000 digit integer

## Problem 9: "Special Pythagorean triplet" Answer: 31875000
def findPythagoreanTripleGivenSum(sum): # returns an array [a,b,c] where a+b+c = sum and a^2+b^2=c^2.
	# without loss of generality, we can say that 1 < a < b < c < sum/2
	# with some algebra, we see that sum*rt(2)/4 < c < sum/2
	# using the bound for c, we see that a < c*rt(2)/2 and b > c*rt(2)/2. Then we have b > sum/4 and a < sum*rt(2)/4
	# code here
# An easy way (paper/pencil) to do this for sum=1000 is to find small Pythagorean triples that sum to a factor of the sum
# and scale each number in the triple by sum/factor. We see that the third Pythagorean triple (for which gcd(a,b)=1) is
# [8,15,17], which sum to 40, which is a factor of 1000. We scale each a, b, and c by 25 and multiply these together.

## Problem 10: "Summation of primes" Answer: 142913828922
def findSumOfPrimes(n): # returns the sum of all primes below n
    import math
	sum = 2
	for i in range(3, n+1, 2):
		sum += i
		for j in range(3, math.floor(i**0.5)+2,2):
			if i % j == 0:
				sum -= i
				break
	return sum

## Problem 11: "Largest product in a grid" Answer: 70600674
def maxProductInGrid(n): # finds maximum product in a p by q grid of n consecutive positive integers
	file = open('general text file.txt') # contains 20x20 grid of integers
	array = []
	for line in file:
		array.append(line.rsplit(" ")) # creates 2D array
	for i in range(0,len(array)): # converts to int
	    for j in range(0,len(array[i])):
			array[i][j].replace("\n", "") # removes new line characters
			array[i][j] = int(array[i][j]) # converts to integers to allow multiplication
	maxVal = -1 # assumes all non-negative integers
	for i in range(0,len(array)):
	    for j in range(0,len(array[i])):
            leftSide = j <= len(array[i])-n
            rightSide = i >= n-1
	        topSide = i <= len(array)-n
	        if topSide: # safe to check down
				product = 1
				for k in range(0, n):
					product *= array[i+k][j]
	            maxVal = max(maxVal, product)
	        if leftSide: # safe to check right
				product = 1
				for k in range(0, n):
					product *= array[i][j+k]
	            maxVal = max(maxVal, product)
	        if topSide and leftSide: # safe to check diagonally: \
				product = 1
				for k in range(0, n):
					product *= array[i+k][j+k]
	            maxVal = max(maxVal, product)
	        if topSide and rightSide: # safe to check diagonally: /
				product = 1
				for k in range(0, n):
					product *= array[i+k][j-k]
	            maxVal = max(maxVal, product)
	print(maxVal)

## Problem 12: "Highly divisible triangular number" Answer: 76576500
def generatePrimes(n, primes=[2, 3]):
    if n < 2:
        return []
    if n == 2:
        return [2]
    for i in range(primes[len(primes)-1]+2, n+1, 2):
        isPrime = True
        for j in primes:
            if i % j == 0:
                isPrime = False
                break
        if isPrime:
            primes.append(i)
    return primes
def factor(n, primes=[]):
    if primes == []:
        primes = generatePrimes(100) # can change default starting point
    factors = {}
    beginAt = 0
    while n != 1:
        for i in range(beginAt, len(primes)):
            check = primes[i]
            if check > n:
                return factors
            while n % check == 0:
                if check in factors:
                    factors[check] += 1
                else:
                    factors[check] = 1
                n /= check
        beginAt = len(primes)
        primes = generatePrimes(min(primes[len(primes)-1] * 2, n), primes)
    return factors
def numFactors(n):
    factors = factor(n)
    product = 1
    for i in factors:
        product *= factors[i]+1
    return product
for i in range(2, 100000): # assumes answer is less than 100000 (could be replaced with while loop)
    if i % 2 == 0:
        value = numFactors(i/2) * numFactors(i+1)
    else:
        value = numFactors(i) * numFactors((i+1)/2)
    if value > 500:
        print(int(i*(i+1)/2))
        break
## Problem 13: "Large sum" Answer: 5537376230
def findFirstDigits(n): # returns the first n digits of the sum of 100 large numbers
	file = open("General text file.txt", "r")
	sum = 0
	for line in file:
		sum += int(line[:(n+3)])
	file.close()
	return int(str(sum)[:n])

## Problem 14: "Longest Collatz sequence" Answer: 837799
def collatzSeqLength(n): # returns the length of the collatz sequence of n
	if n == 1: # base case of recursion
		return 1
	if n%2 == 0: # even case simplification
		return collatzSeqLength(n/2)+1
	return collatzSeqLength((3*n+1)/2)+2 # odd case simplification
def longestCollatzSeq(n): # returns the integer under n generating the longest Collatz sequence (naive implementation)
	maxLength = 1
	maxNum = 1
	for i in range(2,n):
		length = collatzSeqLength(i)
		if length > maxLength:
			maxLength = length
			maxNum = i
	return maxNum

## Problem 15: "Lattice paths" Answer: 137846528820
def combination(n, m): # returns ways to choose m objects out of n objects, independent of order
	import math
	return math.factorial(n)/(math.factorial(m) * math.factorial(n-m))
combination(40, 20)

## Problem 18: finds maximal path in triangular text
def findMaxPath():
	file = open('general text file.txt')
	array = []
	for line in file:
		array.append(line.rsplit(" "))
	array.reverse()
	while len(array) > 1:
			for i in range(0, len(array[1])):
					array[1][i] = int(array[1][i]) + max(int(array[0][i]), int(array[0][i+1]))
			array.pop(0)
	return array[0][0]
  
# Problem 19: "Counting Sundays" Answer: 171
def days_in_year(year):
	if year%4 == 0 and not (year%100 == 0 and year%400 != 0): # leap year
		return 366
	return 365 # not leap year
def days_in_months(year):
	return [31, days_in_year(year)-337, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
def first_of_year(year):
	assert(year >= 1900)
	day_count = 1
	for i in range(1900,year):
		day_count += days_in_year(i)%7
	return day_count%7
def day_of_year(day, month, year, month_length=None):
	day_count = day
	if month_length is None:
		month_length = days_in_months(year)
    for i in range(1, month):
		day_count += month_length[i-1]
	return day_count
def day_of_week(day, month, year, month_length=None):
	assert(year >= 1900)
	if month_length is None:
		month_length = days_in_months(year)
	first_day = first_of_year(year)
	first_day += day_of_year(day, month, year, month_length)-1
	return first_day%7
def prob19():
	sunday_count = 0
	for year in range(1901, 2001):
		month_length = days_in_months(year)
		for month in range(1, 13):
			if day_of_week(1, month, year, month_length) == 0:
				sunday_count += 1
	return sunday_count

## Problem 20: "Factorial digit sum" Answer: 648
def sumOfDigits(n): # returns sum of digits of a number n by converting it to a string and adding each character
	n = str(n)
	sumOfNum = 0
	for i in range(0, len(n)):
		sumOfNum += int(n[i])
	return sumOfNum
import math; n = math.factorial(100); sumOfDigits(n)

## Problem 21: "Amicable numbers" Answer: 31626
def sumOfFactors(n): # returns sum of factors of n
    factorSum = 1
    primeFactorization = factor(n)
    for prime in primeFactorization:
        factorSum *= (prime**(primeFactorization[prime]+1)-1)/(prime-1)
    return int(factorSum)
def amicableNumberSum(n):
    amicableSum = 0
    for i in range(2,n):
        a = sumOfFactors(i)-i
        if i != a and a < n and i == sumOfFactors(a)-a:
            amicableSum += a + i
    return int(amicableSum/2)

## Problem 22: "Names scores" Answer: 871198282
file = open('general text file.txt')
names = []
for line in file:
    names.append(line.rsplit('","'))
names[0][0] = names[0][0][1:]
names[0][len(names[0])-1] = names[0][len(names[0])-1][:-2]
names = names[0]
names.sort()
ans = 0
for i in range(0, len(names)):
    stringSum = 0
    for j in range(0, len(names[i])):
        stringSum += ord(names[i][j])-ord('A')+1
    ans += stringSum * (i+1)
print(ans)

## Problem 23: "Non-abundant sums" Answer: 4179871
def abundantList(n):
    abundantSet = set()
    for i in range(2,n+1):
        if sumOfFactors(i) > 2*i:
            abundantSet.add(i)
    return abundantSet
def sumOf(n, array):
    for num in array:
        if n - num < 0:
            return False
        if n - num in array:
            return True
    return False
def prob23(n):
    nonAbundantSum = 0
    abundantSet = abundantList(n)
    for i in range(1,n+1):
        if not sumOf(i, abundantSet):
            nonAbundantSum += i
    return nonAbundantSum

## Problem 26: "Reciprocal cycles" Answer: 983
def cycle(n):
    num = 2
    maxCycle = 1
    for i in range(2, n):
        if i % 2 != 0 and i % 5 != 0:
            remainders = [1]
            newNum = 0
            isFinished = False
            while not isFinished:
                newNum = remainders[len(remainders)-1] * 10 % i
                isFinished = newNum in remainders
                remainders.append(newNum)
            if len(remainders) > maxCycle:
                num = i
                maxCycle = len(remainders)
    return num

## Problem 27: "Quadratic primes" Answer: -59231
def isPrime(n):
    if n == 2:
        return True
    if n <= 1 or n % 2 == 0:
        return False
    for i in range(3,int(n**0.5)+1,2):
        if n % i == 0:
            return False
    return True
def quadraticPrimes(n): # prints best quadratic prime generator with |a|, |b| < n
    primeList = generatePrimes(n)
    max_count = 0
    best_product = 0
    for a in range(1-n,n,2):
        for b in primeList:
            if b >= n:
                break
            i = 0
            while i**2+a*i+b > 0 and isPrime(i**2+a*i+b):
                i += 1
            if i > max_count:
                max_count = i
                best_product = a*b
    return best_product
	
## Problem 29: "Distinct powers" Answer: 9183
def distinctPowers(n, m): # boring solution
    uniqueSet = set()
    for i in range(n, m+1):
        for j in range(n, m+1):
            uniqueSet.add(i**j)
    return len(uniqueSet)
	
## Problem 30: "Digit fifth powers" Answer: 443839
def digitFifthPowers(n,m): # run with n = 2, m = 300000 (can be reduced to 300000)
    totalSum = 0
    for i in range(n,m):
        currentSum = 0
        temp = i
        while temp > 0:
            currentSum += (temp%10)**5
            temp = floor(temp/10)
        if currentSum == i:
            totalSum += i
    return totalSum
	
## Problem 32: "Pandigital products" Answer: 45228
def pandigital(): # can't be generalized?
	validCombos = set([])
	# start by noticing that we have two possibilites: 1-digit * 4-digit = 4-digit or 2-digit by 3-digit = 4-digit
	# 1-digit case: cannot have first digit be 1 or 9. Simple casework shows cannot be 8
	for i in range(2,8):
		for j in range(1234,9877):
			if hasUniqueDigits([i,j,i*j]):
				validCombos.add(i*j)
	for i in range(12, 98):
		for j in range(123,988):
			if hasUniqueDigits([i,j,i*j]):
				validCombos.add(i*j)
	sum = 0
	for i in validCombos:
		sum += i
	return sum
def hasUniqueDigits(array): # checks if numbers in array are duplicated
	digitsUsed = [0]
	for i in array:
		while i > 0:
			digitsUsed.append(i%10)
			i = (i - i%10)/10
	return len(digitsUsed) == len(set(digitsUsed))

## Problem 34: "Digit factorials" Answer: 40730
def factorialDigitSum(n):
    import math
    allNumbersSum = 0
    for i in range(3,n):
        if i == factorialSum(i):
            allNumbersSum += i
    return allNumbersSum
def factorialSum(n):
    totalSum = 0
    while n > 0:
        totalSum += math.factorial(n%10)
        n = floor(n/10)
    return totalSum

## Problem 44: "Pentagon numbers" Answer: 5482660
def checkPentagon(n):
	a = int((1+(1+24*n)**0.5)/6)
	if (3*a**2-a)/2 == n:
		return True
	return False
def pentagonPair(maxCheck):
	n = 1
	while n < maxCheck:
		for m in range(1,n):
			a = pentagon(n)
			b = pentagon(m)
			if checkPentagon2(a-b) and checkPentagon2(a+b):
				return a-b
		n += 1
	return "Cannot find a pentagonal pair for n,m < " + str(maxCheck)
pentagonPair(10**4)

## Problem 45: "Triangular, pentagonal, and hexagonal" Answer: 1533776805
for i in range(144, 10000000):
	hexagonal = i * (2 * i - 1)
	sol = quadratic(3, -1, -2*hexagonal)
	sol1 = math.floor(sol)
	sol2 = math.ceil(sol)
	if sol1 * (3 * sol1 - 1) == 2*hexagonal or sol2 * (3 * sol2 - 1) == 2*hexagonal:
		print(hexagonal)
		break

## Problem 48: "Self powers" Answer: 9110846700
def selfPowerSum(n):
	count = 0
	for i in range(1,n+1):
		count += i**i
	return count%(10**10)

## Problem 49: "Prime permutations" Answer: 296962999629
def isPermutation(a,b):
    return digitsInNum(a) == digitsInNum(b)
def digitsInNum(n):
    digits = {}
    for i in range(0,10):
        digits[i] = 0
    while n > 0:
        digits[n%10] += 1
        n = floor(n/10)
    return digits
def arithmeticSeqPrimes(n, m):
    primeList = generatePrimes(m)
    primeSet = set(primeList)
    possibleList = []
    for i in primeList:
        if i >= n and i < m:
            for j in primeList:
                if j > i and j < m and isPermutation(i, j):
                    middleNum = int((i+j)/2)
                    if isPermutation(i, middleNum) and middleNum in primeSet:
                        print(i, middleNum, j)

## Problem 53: "Combinatoric selections" Answer:
import math
def combination(n,m): # returns number of ways to choose m numbers from n options
	return math.factorial(n)/(math.factorial(m) * math.factorial(n-m))
def numOfCombinations(n, maximum):
	for i in range(1, int(n/2)+1):
		if combination(n,i) > maximum:
			return n-2*i+1
	return 0
def numOfCombinationsRange(n, maximum):
	sum = 0
	for i in range(1, maximum+1):
		sum += numOfCombinations(i, maximum)
	return sum

## Problem 89: "Roman numerals" Answer: 743
def roman(num):
	if num == 0:
		return ""
	if num >= 1000:
		return "M"+roman(num-1000)
	if num >= 900:
		return "CM"+roman(num-900)
	if num >= 500:
		return "D"+roman(num-500)
	if num >= 400:
		return "CD"+roman(num-400)
	if num >= 100:
		return "C"+roman(num-100)
	if num >= 90:
		return "XC"+roman(num-90)
	if num >= 50:
		return "L"+roman(num-50)
	if num >= 40:
		return "XL"+roman(num-40)
	if num >= 10:
		return "X"+roman(num-10)
	if num >= 9:
		return "IX"+roman(num-9)
	if num >= 5:
		return "V"+roman(num-5)
	if num >= 4:
		return "IV"+roman(num-4)
	if num > 0:
		return "I"+roman(num-1)
def romanToNum(roman):
	if roman[:1] == "M":
		return 1000 + romanToNum(roman[1:])
	if roman[:1] == "D":
		return 500 + romanToNum(roman[1:])
	if roman[:1] == "L":
		return 50 + romanToNum(roman[1:])
	if roman[:1] == "V":
		return 5 + romanToNum(roman[1:])
	if len(roman) != 1:
		if roman[:2] == "CM":
			return 900 + romanToNum(roman[2:])
		if roman[:2] == "CD":
			return 400 + romanToNum(roman[2:])
		if roman[:2] == "XC":
			return 90 + romanToNum(roman[2:])
		if roman[:2] == "XL":
			return 40 + romanToNum(roman[2:])
		if roman[:2] == "IX":
			return 9 + romanToNum(roman[2:])
		if roman[:2] == "IV":
			return 4 + romanToNum(roman[2:])
	if roman[:1] == "C":
		return 100 + romanToNum(roman[1:])
	if roman[:1] == "X":
		return 10 + romanToNum(roman[1:])
	if roman[:1] == "I":
		return 1 + romanToNum(roman[1:])
	return 0
def readRomanNumFile():
	file = open("general text file.txt", "r")
	length = 0
	for line in file:
		line = line[:len(line)-1]
		length += len(line)
		length -= len(roman(int(romanToNum(str(line)))))
	file.close()
	return length

## Problem 92: "Square digit chains" Answer:
def storeDigits(n): # stores digits of a number in an array
	array = []
	while n > 0:
		addVal = n % 10
		if addVal != 0:
			array.append(addVal)
		n = (n - addVal%10)/10
	return array
def digitSquareSum(n):
	count = 0
	for i in range(1, n):
		if squareDigitChainEnd(i) == 89:
			count += 1
	return count
def squareDigitChainEnd(n):
	array = storeDigits(n)
	sum = 0
	for j in range(0, len(array)):
		sum += array[j]**2
	if sum == 1:
		return 1
	if sum == 89:
		return 89
	return squareDigitChainEnd(sum)

## Problem 97: "Large non-Mersenne prime" Answer: 8739992577
def prob97(n):
	a = n % 40
	div = (n - n%40)/40
	mult = (28433 * 2**a) % 10**10
	b = 2**40 % 10**10
	for i in range(0, int(div)):
		mult *= b
		mult = mult % 10**10
	return mult+1
	# a takes the modularity of n
###########
def primeDays():
    month_length = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    for i in range(1,13):
        for j in range(1,10):
            if isPrime(int(str(i)+"0"+str(j)+"2017")): print(str(i)+"/0"+str(j)+"/2017")
        for j in range(10,month_length[i-1]+1):
            if isPrime(int(str(i)+str(j)+"2017")): print(str(i)+"/"+str(j)+"/2017")