import asyncio
import os

import numpy as np
from PIL import Image
from playwright.async_api import async_playwright

os.makedirs('src/ui/vtk_viewer/tests/artifacts', exist_ok=True)
SWIFTSHADER = ["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader"]

def non_black(path):
    rgb = np.asarray(Image.open(path).convert("RGB"), dtype=np.int32)
    return float((rgb.sum(axis=2) > 15).mean())

async def main():
    errs = []
    async with async_playwright() as p:
        b = await p.chromium.launch(headless=True, args=SWIFTSHADER)
        page = await b.new_page(viewport={'width': 1200, 'height': 900})
        page.on("console", lambda m: errs.append(f"{m.type}: {m.text}") if m.type in ('error', 'warning') else None)
        page.on("pageerror", lambda e: errs.append("PAGEERR: " + str(e)))
        await page.goto("http://localhost:8501", wait_until="domcontentloaded", timeout=60000)
        await page.wait_for_timeout(5000)
        vf = None
        for f in page.frames:
            if 'srcdoc' in f.url:
                vf = f
                break
        print("srcdoc viewer frame:", vf is not None)
        if vf:
            ok = False
            for _ in range(240):
                try:
                    st = await vf.evaluate("""() => {
                        const c=document.querySelector('canvas');
                        return {
                            c:c?{w:c.width,h:c.height}:null,
                            msg:document.getElementById('msg').textContent,
                            load:document.getElementById('loading').style.display,
                            controls:document.getElementById('controls').childElementCount,
                        };
                    }""")
                except Exception as e:
                    st = {"err": str(e)}
                if st.get('c') and st['c']['w'] > 0 and st.get('load') in ('none', 'hidden', ''):
                    ok = True
                    break
                if st.get('msg'):
                    print("MSG:", st['msg'][:300])
                    break
                await page.wait_for_timeout(500)
            print("render ok:", ok, st if not ok else "")
            if ok:
                await page.wait_for_timeout(2000)
                diag = await vf.evaluate("window.__viewer_diag(false)")
                print('served runtime:', {k: diag.get(k) for k in (
                    'implementation', 'volumeCount', 'componentCount',
                    'independentComponents', 'autoAdjust', 'sampleDistance',
                    'shade', 'dpr', 'domSize', 'cssSize')})
                p2d = 'src/ui/vtk_viewer/tests/artifacts/demo_served_2d.png'
                p3d = 'src/ui/vtk_viewer/tests/artifacts/demo_served_3d.png'
                await vf.locator('#slice-canvas').screenshot(path=p2d)
                await vf.locator('#viewer').screenshot(path=p3d)
                print('served split non-black:', {'2d': non_black(p2d), '3d': non_black(p3d)})
                box = await page.main_frame.evaluate("""() => {
                    const els=document.getElementsByTagName('iframe');
                    for(const e of els){ const r=e.getBoundingClientRect(); if(r.width>50) return {x:r.x,y:r.y,w:r.width,h:r.height}; }
                    return null;
                }""")
                print("iframe box:", box)
                if box:
                    clip = {"x": box['x'], "y": box['y'], "width": min(box['w'], 1200), "height": min(box['h'], 850)}
                    await page.screenshot(path="src/ui/vtk_viewer/tests/artifacts/volume_single_channel.png", clip=clip)
                await page.screenshot(path="src/ui/vtk_viewer/tests/artifacts/page.png")
        print("frames:", [f.url for f in page.frames])
        print("=== errs ===")
        for e in errs:
            if 'Unrecognized' in e or 'same-origin' in e:
                continue
            print(e)
        await b.close()

asyncio.run(main())
