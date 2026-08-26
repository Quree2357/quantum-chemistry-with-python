# 8.2. 지름 방정식과 양자수

앞 절에서 수소 원자의 Schrödinger 방정식을 살펴보았습니다. 그리고 변수분리를 통해서 지름 성분과 각도 성분의 방정식을 따로 얻었죠. 각도 성분의 파동 방정식의 해는 7장에서 봤던 구면 조화 함수로 나타났습니다. 이제 지름 성분의 파동 방정식을 풀 차례입니다. 복잡하게 생겼지만 천천히 하나씩 해결해나가면 풀 수 있습니다.


## 지름 방정식

수소 원자의 지름 방정식을 다시 한번 살펴봅시다.

$$
-\frac{\hbar^2}{2m_e r^2}\frac{d}{dr} \left( r^2 \frac{d}{dr} \right) R(r) + \left[ \frac{l(l+1)\hbar^2}{2m_e r^2} - \frac{e^2}{4\pi \varepsilon_0 r} \right] R(r) = ER(r)
$$

다시 봐도 끔찍하게 생겼군요... 일단 맨 왼쪽 항에 있는 이중 미분부터 처리해봅시다. 이렇게 치환을 해보면 어떨까요? $u(r) = rR(r)$. 한번 대입해봅시다.

$$
\frac{1}{r^2}\frac{d}{dr}\left(r^2\frac{dR}{dr}\right)=\frac{1}{r^2}\frac{d}{dr}\left(r^2\frac{d}{dr}\frac{u}{r}\right)=\frac{1}{r^2}\frac{d}{dr}\left(r^2 \left(\frac{1}{r}\frac{du}{dr}-\frac{u}{r^2}\right)\right)=\frac{1}{r^2}\frac{d}{dr}\left(r\frac{du}{dr}-u\right)
=\frac{1}{r^2}\left(r\frac{d^2u}{dr^2}\right)=\frac{1}{r}\frac{d^2u}{dr^2}
$$

치환하고 보니 복잡한 미분항이 간단한 이계 미분항으로 바뀌었습니다! 원래 지름 방정식에 넣고 전체에 $r$을 곱해서 정리하면 다음과 같이 되겠네요.

$$
-\frac{\hbar^2}{2m_e}\frac{d^2}{dr^2}u(r) + \left[ \frac{l(l+1)\hbar^2}{2m_e r^2} - \frac{e^2}{4\pi \varepsilon_0 r} \right]u(r) = Eu(r)
$$

앞에서 계속 봐왔던 1차원 Schrödinger 방정식 형태가 되었습니다. 위 식에서 $u(r)$에 곱해진 대괄호 부분을 통째로 $V_{\text{eff}}(r)$라고 하면 완전히 모양이 똑같아집니다.

$$
-\frac{\hbar^2}{2m_e}\frac{d^2}{dr^2}u(r) + V_{\text{eff}}(r)u(r) = Eu(r)
$$

경계조건도 잊지 말고 확인해야죠. $u(r)=rR(r)$이었으니까 $u(0)=0$이어야 하고, 올바른 파동 함수가 되려면 당연히 무한대의 거리에서는 $u(r) \to 0$이어야 합니다.
