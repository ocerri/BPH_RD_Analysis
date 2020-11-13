cat = {}

class Bauble: pass

cat['high'] = Bauble()
cat['high'].name = 'High'
cat['high'].min_pt = 12
cat['high'].max_pt = 1e4
cat['high'].trg = 'Mu12_IP6'
cat['high'].minIP = 7

cat['mid'] = Bauble()
cat['mid'].name = 'Mid'
cat['mid'].min_pt = 9
cat['mid'].max_pt = 12
cat['mid'].trg = 'Mu9_IP6'
cat['mid'].minIP = 6

cat['low'] = Bauble()
cat['low'].name = 'Low'
cat['low'].min_pt = 7
cat['low'].max_pt = 9
cat['low'].trg = 'Mu7_IP4'
cat['low'].minIP = 4

cat['single'] = Bauble()
cat['single'].name = 'Single'
cat['single'].min_pt = 9
cat['single'].max_pt = 1e6
cat['single'].trg = 'Mu9_IP6'
cat['single'].minIP = 6

categories = cat
