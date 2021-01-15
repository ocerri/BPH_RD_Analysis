cat = {}

class Bauble: pass

cat['high'] = Bauble()
cat['high'].name = 'High'
cat['high'].min_pt = 12.2
cat['high'].max_pt = 1e4
cat['high'].trg = 'Mu12_IP6'
cat['high'].minIP = 7

cat['mid'] = Bauble()
cat['mid'].name = 'Mid'
cat['mid'].min_pt = 9.2
cat['mid'].max_pt = 12.2
cat['mid'].trg = 'Mu9_IP6'
cat['mid'].minIP = 7

cat['low'] = Bauble()
cat['low'].name = 'Low'
cat['low'].min_pt = 7.2
cat['low'].max_pt = 9.2
cat['low'].trg = 'Mu7_IP4'
cat['low'].minIP = 5

cat['single'] = Bauble()
cat['single'].name = 'Single'
cat['single'].min_pt = 9.2
cat['single'].max_pt = 1e6
cat['single'].trg = 'Mu9_IP6'
cat['single'].minIP = 7

cat['none'] = None

categories = cat
