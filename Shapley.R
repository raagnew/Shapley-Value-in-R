rm(list=ls())
shapley.value <- function(v){
#  Bob Agnew (raagnew1@gmail.com, www.raagnew.com)
#  This function computes Shapley value for a cooperative game
#  in characteristic function form.  n is the number of players.
#  v is the characteristic function with dimension rep(2,n).
#  v[k1,...,kn] is the value of a coalition of players with 
#  argument 2 for inclusion, 1 for exclusion.  v[1,...,1] = 0
#  corresponds to the empty coalition.  v[2,...,2] is the value
#  of the grand coalition consisting of all n players.  Output
#  s  is the Shapley value vector indicating the relative
#  bargaining power and expected share of the pie for each player.
#  Uses intToBits function and is limited to 32 players, although
#  the computational limit is less than that since 2^32 > 4 billion.
#  Reference is G. Owen, Game Theory, 3rd Ed, Academic Press, 2001.
n <- length(dim(v))
v <- as.vector(v)
m <- array(as.integer(intToBits(-1+1:2^n)),dim=c(32,2^n))[1:n,]
t1 <- c(0,((t <- colSums(m))-1)[-1])
w <- factorial(t1)*factorial(n-t)/factorial(n)
s <- NULL
for (i in 1:n){
v0 <- v
v0[(1:2^n)[m[i,]==1]] <- v[-2^(i-1)+(1:2^n)[m[i,]==1]]
s <- c(s,sum(w*(v-v0)))}
s}

#  Example 1
#  Example XII.1.8 in Owen.  We have a corporation with four
#  shareholders having 10, 20, 30, and 40 shares respectively with
#  decisions settled by a simple majority (51%) of shares.  This
#  is a simple zero-one game with minimal winning coalitions
#  { 2, 4 }, { 3, 4 }, and { 1, 2, 3 }.  As Owen points out, the Shapley
#  value vector is ( 1/12, 1/4, 1/4, 5/12 ) which contrasts with the
#  "vote vector" ( 1/10, 1/5, 3/10, 2/5 ) so player 1 has less power
#  than his shares reflect, player 4 has more power than his
#  shares reflect, and players 2 and 3 have equal power, despite
#  player 3's extra shares.  We note that there is an efficient
#  generating function algorithm for computing Shapley value for much
#  larger weighted voting games of this sort, including the U.S. electoral
#  college system.  That approach is discussed in Owen and is widely
#  available in online and downloadable calculators.  
v <- array(0,dim=rep(2,4))
v[,2,,2] <- v[,,2,2] <- v[2,2,2,] <- 1
shapley.value(v)

#  Example 2
#  This is another simple zero-one game with one big player 1 and
#  four little players 2, 3, 4, and 5.  The big player can win with
#  any of the little players and the little players can win if they
#  all cooperate and freeze out the big player.  Shapley value
#  indicates that the big player has 60% of the bargaining power while
#  the little players share the remaining 40% equally.
v <- array(0,dim=rep(2,5))
v[2,2,,,] <- v[2,,2,,] <- v[2,,,2,] <- v[2,,,,2] <- v[,2,2,2,2] <- 1
shapley.value(v)

#  Example 3
#  This is a more complex market game where players 1 and 2 are sellers
#  and players 3, 4, and 5 are buyers.  The product is differentiated
#  in terms of seller's production cost and buyer's valuation (for resale
#  or use), but any seller can satisfy any buyer's demand.  Players 1
#  and 2 have supplies (capacities) of 10 and 20 units respectively.
#  Players 3, 4, and 5 have demands (requirements) of 5, 20, and 10
#  units respectively, so there are 5 units of excess demand.
#  Players 1 and 2 have the following unit production costs to supply
#  players 3, 4, and 5.
#
#       3     4     5 
#  1 | $50 | $50 | $40 |
#  2 | $30 | $20 | $15 |
#
#  Players 3, 4, and 5 have the following unit valuations for supplies
#  from players 1 and 2.
#
#        3     4     5
#  1 | $100 | $95 | $90 |
#  2 |  $80 | $75 | $75 |
#
#  Hence, total unit "spreads" or "margins" between buyers' valuations
#  and sellers' costs are as follows.
#
#       3     4     5
#  1 | $50 | $45 | $50 |
#  2 | $50 | $55 | $60 |
#
#  If we ignore antitrust laws and uniform pricing laws, the joint optimal
#  solution, with all players participating, yields total margin of
#  $1625 to be split among the five players per these unit allocations.
#
#      3    4    5
#  1 | 5 |  5 |  0 |
#  2 | 0 | 10 | 10 |
#
#  However, this allocation, excluding player 3, yields a close $1600.
#
#       4    5
#  1 | 10 |  0 |
#  2 | 10 | 10 |
#
#  Solving all conceivable subproblems for positive margin coaltions
#  yields the following characteristic function values
#
#        {1,3} |  $250
#        {1,4} |  $450
#        {1,5} |  $500
#        {2,3} |  $250
#        {2,4} | $1100
#        {2,5} |  $600
#      {1,2,3} |  $250
#      {1,2,4} | $1100
#      {1,2,5} |  $600
#      {1,3,4} |  $475
#      {1,3,5} |  $500
#      {1,4,5} |  $500
#      {2,3,4} | $1100
#      {2,3,5} |  $850
#      {2,4,5} | $1150
#    {1,2,3,4} | $1350
#    {1,2,3,5} |  $850
#    {1,2,4,5} | $1600
#    {1,3,4,5} |  $500
#    {2,3,4,5} | $1150
#  {1,2,3,4,5} | $1625
#  Assuming contracting and pricing flexibility, Shapley value
#  ( $264.17, $624.58, $72.50, $443.33, $220.42 ) indicates a reasonable
#  allocation of total margin over supplier production cost to the five
#  players.  Obviously, the low-cost supplier (player 2) gets a big share
#  while the limited buyer (player 3) gets a small share.  We can gain
#  additional insight by dissecting this result in terms of unit margins.
#  Player 3 gets $14.50 margin for each of the 5 units he buys from
#  player 1 so that player 1 gets $50 - $14.50  =  $35.50 for each of
#  these units.  Put another way, player 1 charges player 3 $50 + $35.50
#  = $85.50 for each of these units, leaving player 3 with $14.50
#  markup to his $100 valuation.  Similarly, player 5 gets $22.04 margin
#  for each of the 10 units he buys from player 2 so that player 2 gets
#  $60 - $22.04 = $37.96 for each of these units.  Put another way,
#  player 2 charges player 5  $15 + $37.96  =  $52.96 for each of these
#  units, leaving player 5 with $22.04 markup to his $75 valuation.
#  Player 1 gets ($264.17 - 5 * $35.50) / 5  =  $17.33 margin for each
#  of the 5 units he supplies to player 4.  Put another way, he charges
#  player 4 $50 + $17.33  =  $67.33 for each of these units, leaving         
#  player 4 $27.67 markup to his $95 valuation.  Finally, player 2
#  gets ($624.58 - 10 * 37.96) / 10  =  $24.50 for each of the 10 units
#  he supplies to player 4.  Put another way, he charges player 4
#  $20 + $24.50  =  $44.50 for each of these units, leaving player 4
#  with $30.50 markup to his $75 valuation.  This is of course an
#  idealized situation, but Shapley value shows how "fair" margins and
#  prices can be conceptualized in this setting. 
v <- array(0,dim=rep(2,5))
v[2,1,2,1,1] <- 250
v[2,1,1,2,1] <- 450
v[2,1,1,1,2] <- 500
v[1,2,2,1,1] <- 250
v[1,2,1,2,1] <- 1100
v[1,2,1,1,2] <- 600
v[2,2,2,1,1] <- 250
v[2,2,1,2,1] <- 1100
v[2,2,1,1,2] <- 600
v[2,1,2,2,1] <- 475
v[2,1,2,1,2] <- 500
v[2,1,1,2,2] <- 500
v[1,2,2,2,1] <- 1100
v[1,2,2,1,2] <- 850
v[1,2,1,2,2] <- 1150
v[2,2,2,2,1] <- 1350
v[2,2,2,1,2] <- 850
v[2,2,1,2,2] <- 1600
v[2,1,2,2,2] <- 500
v[1,2,2,2,2] <- 1150
v[2,2,2,2,2] <- 1625
shapley.value(v)

#  Example 4
#  This is an international trade game with four countries negotiating a
#  multilateral agreement on behalf of their citizens.  The issue is how
#  the total gain from trade should be divided among the four players with
#  the following GDPs in $ billion.
#
#  1 | Large Developed   | $15000
#  2 | Medium Developed  |  $2000
#  3 | Large Developing  |  $1000
#  4 | Medium Developing |   $400
#
#  Gains from trade in $ billion from various coalitions are as follows.
#
#            1     2    3    4    Total
#      {1,2} |  $60 |  $60 |  $0 |  $0 || $120
#      {1,3} |  $45 |   $0 | $40 |  $0 ||  $85
#      {1,4} |  $15 |   $0 |  $0 | $40 ||  $55
#      {2,3} |   $0 |  $20 | $20 |  $0 ||  $40
#      {2,4} |   $0 |  $12 |  $0 | $16 ||  $28
#      {3,4} |   $0 |   $0 | $10 |  $8 ||  $18
#    {1,2,3} |  $90 |  $80 | $70 |  $0 || $240
#    {1,2,4} |  $75 |  $70 |  $0 | $60 || $205
#    {1,3,4} |  $60 |   $0 | $50 | $52 || $162
#    {2,3,4} |   $0 |  $30 | $25 | $20 ||  $75
#  {1,2,3,4} | $105 | $100 | $80 | $80 || $365  Free trade split
#
#  This last grand coalition split represents the "free trade" solution.
#  But the Shapley "fair trade" solution ( $137.58, $96.58, $74.75, $56.08 )
#  accounts for all the subcoalitions that could form and allocates more of
#  the total $365 gain to the large developed country and less to the developing
#  countries, particularly to the smallest one, via cash transfers, not tariffs
#  on specific products.  Country 1 in effect charges differentiated annual fees
#  for access to its market.  The developing countries still achieve large
#  percentage GDP gains from the Shapley value fair trade solution.
gdp <- c(15000,2000,1000,400)
v <- array(0,dim=rep(2,4))
v[2,2,1,1] <- 120
v[2,1,2,1] <- 85
v[2,1,1,2] <- 55
v[1,2,2,1] <- 40
v[1,2,1,2] <- 28
v[1,1,2,2] <- 18
v[2,2,2,1] <- 240
v[2,2,1,2] <- 205
v[2,1,2,2] <- 162
v[1,2,2,2] <- 75
v[2,2,2,2] <- 365
s <- shapley.value(v)
#  Fair trade split
s
#  Percent of total
100*s/sum(s)
#  Percent GDP gain
100*s/gdp

#  Example 5
#  This a larger version of example 2 which allows us to demonstrate
#  automated input and also to evaluate run time.  There are three
#  big players (1,2,3) and twenty little players (4,...,23).  Any two
#  big players can partner with any two little players to win and all
#  twenty little players can partner to win.  Shapley value indicates
#  that the three big players equally share 95.4% of the bargaining
#  power while the little players share the remaining 4.6% equally.
#  This code runs quickly on an ordinary PC even though there are
#  2^23  =  8,388,608 coalitions.
date()
v <- array(0,dim=rep(2,23))
for (i1 in 1:2){
for (i2 in 1:2){
for (i3 in 1:2){
for (i4 in 1:2){
for (i5 in 1:2){
for (i6 in 1:2){
for (i7 in 1:2){
for (i8 in 1:2){
for (i9 in 1:2){
for (i10 in 1:2){
for (i11 in 1:2){
for (i12 in 1:2){
for (i13 in 1:2){
for (i14 in 1:2){
for (i15 in 1:2){
for (i16 in 1:2){
for (i17 in 1:2){
for (i18 in 1:2){
for (i19 in 1:2){
for (i20 in 1:2){
for (i21 in 1:2){
for (i22 in 1:2){
for (i23 in 1:2){
a <- c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23) - 1
v[i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23] <- (sum(a[1:3]) >= 2)&(sum(a[4:23]) >= 2)
}}}}}}}}}}}}}}}}}}}}}}}
v[,,,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2] <- 1
s <- shapley.value(v)
s
sum(s[1:3])
sum(s[4:23])
date()





